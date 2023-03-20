/*
 * Copyright 2020 Improbable Worlds Limited
 * Copyright 2022 Tilmann Zäschke
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#ifndef PHTREE_V20_NODE_H
#define PHTREE_V20_NODE_H

#include "entry.h"
#include "phtree/common/b_plus_tree_hash_map.h"
#include "phtree/common/b_plus_tree_multimap.h"
#include "phtree/common/common.h"
#include <iostream>
#include <map>

namespace improbable::phtree::v20 {

/*
 * We provide different implementations of the node's internal entry set.
 * All implementations are equivalent to "std::map<hc_pos_t, Entry>" which can be used as
 * a plugin example for verification.
 *
 * - `array_map` is the fastest, but has O(2^DIM) space complexity. This can be very wasteful
 *   because many nodes may have only 2 entries.
 *   Also, iteration depends on some bit operations and is also O(DIM) per step if the CPU/compiler
 *   does not support CTZ (count trailing zeroes).
 * - `sparse_map` is slower, but requires only O(n) memory (n = number of entries/children).
 *   However, insertion/deletion is O(n), i.e. O(2^DIM) time complexity in the worst case.
 * - 'b_plus_tree_map` is the least efficient for small node sizes but scales best with larger
 *   nodes and dimensionality. Remember that n_max = 2^DIM.
 *
 *   For the LN (large-node) implementation, we require the following logic:
 *   - emplace(): succeeds, unless:
 *     a) There is a subnode with the same hc_pos & prefix
 *     b) There is another entry with the same hc_pos & key
 *   - lower_bound(hc_pos):
 *     a) return iterator to first entry with same hc_pos
 *   - lower_bound(hc_pos, key):
 *     a) return iterator to first entry with same hc_pos && subnode/prefix
 *     b) return iterator to first entry with same hc_pos && key
 *   - erase(hc_pos, key)
 *     a) erase entry with matching hc_pos & key
 *     [NOT REQUIRED: b) erase entry with matching hc_pos & subnode/prefix]
 *
 *   We can avoid the key/prefix matches by manually iterating over entries (which is already
 * implemented). So we only need:
 *   - lower_bound(key)
 *   - emplace_hint(hint, key, value)
 *   - erase(hint)
 *
 *
 *   TODO
 *   MM rules:
 *   - A quadrants contains either a subnode OR value-entries, never both.
 *   - "over-flow" entries are added until node-size becomes MAX. After that, new entries are
 *     only added if their quadrant is empty. If the quadrant contains an entry, a subnode is
 * created. if it contains a subnode, the subnode is traversed (as always when subnodes are
 * encountered).
 *     --> This speeds up look-up and sub-node navigation. Multiple entries are only traversed in
 * leaf-quadrants.
 *
 *  New concept: Leaf quadrants vs inner (=branch?) quadrants.
 *
 *
 */
template <dimension_t DIM, typename Entry>
using EntryMap =
    typename std::conditional_t < DIM<32,
                                      ::phtree::bptree::b_plus_tree_multimap<std::uint32_t, Entry>,
                                      ::phtree::bptree::b_plus_tree_multimap<std::uint64_t, Entry>>;

template <dimension_t DIM, typename Entry>
using EntryIterator = typename std::remove_const_t<decltype(EntryMap<DIM, Entry>().begin())>;
template <dimension_t DIM, typename Entry>
using EntryIteratorC = decltype(EntryMap<DIM, Entry>().cbegin());

/*
 * A node of the PH-Tree. It contains up to 2^DIM entries, each entry being either a leaf with data
 * of type T or a child node (both are of the variant type Entry).
 *
 * The keys (coordinates) of all entries of a node have the same prefix, where prefix refers to the
 * first 'n' bits of their keys. 'n' is equivalent to "n = w - GetPostLen() - 1", where 'w' is the
 * number of bits of the keys per dimension (usually w = 64 for `int64_t` or 'double').
 *
 * The entries are stored in an EntryMap indexed and ordered by their "hypercube address".
 * The hypercube address is the ID of the quadrant in the node. Nodes are effectively binary
 * hypercubes (= binary Hamming space) on {0,1}^DIM. The hypercube address thus uses one bit per
 * dimension to address all quadrants of the node's binary hypercube. Each bit designates for one
 * dimension which quadrant it refers to, such as 0=left/1=right; 0=down/1=up; 0=front/1=back; ... .
 * The ordering of the quadrants thus represents a z-order curve (please note that this completely
 * unrelated to `z-ordering` used in graphics).
 *
 * A node always has at least two entries, except for the root node which can have fewer entries.
 * None of the functions in this class are recursive, see Emplace().
 */
template <dimension_t DIM, typename T, typename SCALAR>
class Node {
    using bit_width_t = detail::bit_width_t;
    using KeyT = PhPoint<DIM, SCALAR>;
    using EntryT = Entry<DIM, T, SCALAR>;
    using hc_pos_t = detail::hc_pos_dim_t<DIM>;

    static_assert(std::is_move_assignable_v<EntryMap<DIM, EntryT>>);
    static_assert(std::is_move_constructible_v<EntryMap<DIM, EntryT>>);

    // static constexpr hc_pos_t MAX_SIZE = std::max(hc_pos_t(8), hc_pos_t(1) << DIM);
    // TODO test with 2D, use something like std::max(f(DIM), 8) ?
    // TODO think about limiting this for high DIM, nodes with 2^DIM entries are not good :-D
    static constexpr hc_pos_t NODE_MAX = (hc_pos_t(1) << DIM) + hc_pos_t(0);

  public:
    Node() : entries_{} {}

    // Nodes should (almost) never be copied!
    Node(const Node&) = delete;
    Node(Node&&) noexcept = default;
    Node& operator=(const Node&) noexcept = delete;
    Node& operator=(Node&&) noexcept = default;

    [[nodiscard]] auto GetEntryCount() const {
        return entries_.size();
    }

    /*
     * Attempts to emplace an entry in this node.
     * The behavior is analogous to std::map::emplace(), i.e. if there is already a value with the
     * given hypercube address 'hc_pos', that value is returned. This function is also
     * non-recursive, it will return a child node instead of traversing it.
     *
     * The scenarios in detail:
     *
     * If there is no entry at the position of 'hc_pos', a new entry is inserted. The new entry is
     * constructed from a constructor of T that accepts the arguments `args`. Also, 'is_inserted' is
     * set top 'true'.
     *
     * If there is an entry with a value T at 'hc_pos', that value is returned. The value is
     * _not_ overwritten.
     *
     * If there is a child node at the position of 'hc_pos', the child node's prefix is analysed.
     * If the prefix indicates that the new value would end up inside the child node or any of its
     * children, then the child node is returned for further traversal.
     * If the child nodes' prefix is different, then a new node is created. The new node contains
     * the child node and a new key/value entry constructed with `args`. The new node is inserted in
     * the current node at the position of the former child node. The new value is returned and
     * 'is_inserted' is set to 'true'.
     *
     * @param is_inserted The function will set this to true if a new value was inserted
     * @param key The key for which a new value should be inserted.
     * @param args Constructor arguments for creating a value T that can be inserted for the key.
     */
    template <typename... Args>
    EntryT& Emplace(bool& is_inserted, const KeyT& key, EntryT* parent_entry, Args&&... args) {
        bit_width_t postfix_len = parent_entry->GetNodePostfixLen();
        hc_pos_t hc_pos = detail::CalcPosInArray(key, postfix_len);
        // subnode? return for further traversal
        auto iter = entries_.lower_bound(hc_pos);
        if (postfix_len > 0 && iter != entries_.end() && iter->first == hc_pos) {
            if (iter->second.IsNode()) {
                return HandleCollision(
                    iter->second, is_inserted, key, postfix_len, std::forward<Args>(args)...);
            }

            // TODO DO math for this splitting policy.
            //  Node may contain 8 + 2^(d-1) entries max....

            // split node if it is too large
            if (entries_.size() >= NODE_MAX) {
                return SplitQuadrantAndInsert(
                    iter, is_inserted, postfix_len, key, parent_entry, std::forward<Args>(args)...);
            }
        }

        // Does not exist
        is_inserted = true;
        auto entry_iter =
            entries_.emplace_hint(iter, hc_pos, EntryT{key, std::forward<Args>(args)...});
        return entry_iter->second;
    }

    /*
     * Use this when we know
     * - that there are no duplicates (skip comparing keys)
     * - and that there is going to be no overflow (skip counting quadrant)
     */
    template <typename... Args>
    EntryT& EmplaceUnchecked(
        bool& is_inserted, const KeyT& key, bit_width_t postfix_len, Args&&... args) {
        hc_pos_t hc_pos = detail::CalcPosInArray(key, postfix_len);

        // subnode? return for further traversal
        auto iter = entries_.lower_bound(hc_pos);
        if (iter != entries_.end() && iter->second.IsNode() && iter->first == hc_pos) {
            return HandleCollision(
                iter->second, is_inserted, key, postfix_len, std::forward<Args>(args)...);
        }

        // Does not exist
        is_inserted = true;
        auto entry_iter =
            entries_.emplace_hint(iter, hc_pos, EntryT{key, std::forward<Args>(args)...});
        return entry_iter->second;
    }

    /*
     * Returns the value (T or Node) if the entry exists and matches the key. Child nodes are
     * _not_ traversed.
     * @param key The key of the entry
     * @param parent The parent node
     * @return The sub node or null.
     */
    // TODO is this function still useful?
    EntryT* Find(const KeyT& key, bit_width_t postfix_len) {
        hc_pos_t hc_pos = CalcPosInArray(key, postfix_len);
        auto iter = entries_.find(hc_pos);
        if (iter != entries_.end() && iter->second.IsNode()) {
            if (DoesEntryMatch(iter->second, key, postfix_len)) {
                return &iter->second;
            }
            return nullptr;
        }
        while (iter != entries_.end() && iter->first == hc_pos) {
            if (iter->second.GetKey() == key) {
                return &iter->second;
            }
            ++iter;
        }
        return nullptr;
    }

    const EntryT* FindC(const KeyT& key, bit_width_t postfix_len) const {
        return const_cast<Node&>(*this).Find(key, postfix_len);
    }

    // TODO is this the same as LowerBound?
    auto FindIter(const KeyT& key, bit_width_t postfix_len) {
        hc_pos_t hc_pos = detail::CalcPosInArray(key, postfix_len);
        // TODO use lower_bound?
        auto iter = entries_.find(hc_pos);
        if (iter != entries_.end() && iter->second.IsNode()) {
            if (DoesEntryMatch(iter->second, key, postfix_len)) {
                return iter;
            }
            return entries_.end();
        }
        while (iter != entries_.end() && iter->first == hc_pos) {
            if (iter->second.GetKey() == key) {
                return iter;
            }
            ++iter;
        }
        return entries_.end();
    }

    // TODO use this in Erase() ?!?!
    auto LowerBound(const KeyT& key, bit_width_t postfix_len, bool& found) {
        hc_pos_t hc_pos = detail::CalcPosInArray(key, postfix_len);
        auto iter = entries_.lower_bound(hc_pos);
        while (iter != entries_.end() && iter->first == hc_pos) {
            if (iter->second.IsNode()) {
                found = DoesEntryMatch(iter->second, key, postfix_len);
                // TODO is this okay? The client cannot tell that this is wrong...??!?
                return iter;
            }
            if (iter->second.GetKey() == key) {
                found = true;
                return iter;
            }
            ++iter;
        }
        found = false;
        return iter;
    }

    auto LowerBound(const KeyT& key, bit_width_t postfix_len, bool& found) const {
        return const_cast<Node&>(*this).LowerBound(key, postfix_len, found);
    }

    auto LowerBound(const KeyT& key, const T& value, bit_width_t postfix_len, bool& found) {
        hc_pos_t hc_pos = detail::CalcPosInArray(key, postfix_len);
        auto iter = entries_.lower_bound(hc_pos);
        if (iter != entries_.end() && iter->second.IsNode()) {
            found = DoesEntryMatch(iter->second, key, postfix_len);
            return iter;
        }
        while (iter != entries_.end() && iter->first == hc_pos) {
            if (iter->second.GetValue() == value && iter->second.GetKey() == key) {
                found = true;
                return iter;
            }
            ++iter;
        }
        found = false;
        return iter;
    }

    auto LowerBound(const KeyT& key, const T& value, bit_width_t postfix_len, bool& found) const {
        return const_cast<Node&>(*this).LowerBound(key, value, postfix_len, found);
    }

    auto End() {
        return entries_.end();
    }

    auto End() const {
        return entries_.end();
    }

    auto FindPrefix(
        const KeyT& prefix, bit_width_t prefix_post_len, bit_width_t node_postfix_len) const {
        assert(prefix_post_len <= node_postfix_len);
        hc_pos_t hc_pos = detail::CalcPosInArray(prefix, node_postfix_len);
        const auto iter = entries_.find(hc_pos);
        if (iter == entries_.end() || iter->second.IsValue() ||
            iter->second.GetNodePostfixLen() < prefix_post_len) {
            // We compare the infix only if it lies fully within the prefix.
            return entries_.end();
        }

        if (DoesEntryMatch(iter->second, prefix, node_postfix_len)) {
            return EntryIteratorC<DIM, EntryT>{iter};
        }
        return entries_.end();
    }

    /*
     * Attempts to erase a key/value pair.
     * This function is not recursive, if the 'key' leads to a child node, the child node
     * is returned and nothing is removed.
     *
     * @param key The key of the key/value pair to be erased
     * @param parent_entry The parent node of the current node (=nullptr) if this is the root node.
     * @param allow_move_into_parent Whether the node can be merged into the parent if only 1
     * entry is left.
     * @param found This is and output parameter and will be set to 'true' if a value was removed.
     * @return A child node if the provided key leads to a child node.
     */
    EntryT* Erase(const KeyT& key, const T& value, EntryT* parent_entry, size_t& found) {
        auto postfix_len = parent_entry->GetNodePostfixLen();
        hc_pos_t hc_pos = detail::CalcPosInArray(key, postfix_len);
        auto it = entries_.find(hc_pos);
        if (it != entries_.end()) {
            if (it->second.IsNode()) {
                if (DoesEntryMatch(it->second, key, postfix_len)) {
                    return &it->second;
                }
                return nullptr;
            }

            // TODO document or fix: we are removing only the FIRST occurrence!!!
            //  Fix: implement std::erase_if() for b_plus_tree, see "Notes" on erase-remove idiom:
            //    https://en.cppreference.com/w/cpp/algorithm/remove
            while (it != entries_.end() && it->first == hc_pos) {
                // TODO can we skip the first key comparison?
                //    -> move DoesEntryMatch from above into if-clause with "if (IsNode)"!
                if (it->second.GetValue() == value && it->second.GetKey() == key) {
                    entries_.erase(it);
                    found = 1;
                    CheckUnderflow(parent_entry);
                    return nullptr;
                }
                ++it;
            }
        }
        return nullptr;
    }

    void Erase(
        EntryIteratorC<DIM, EntryT> iter, EntryT* parent_entry, bool check_underflow = true) {
        assert(iter != entries_.end());
        assert(iter->second.IsValue());
        entries_.erase(iter);
        if (check_underflow) {
            CheckUnderflow(parent_entry);
        }
    }

    void CheckUnderflow(EntryT* parent_entry) {
        if (GetEntryCount() == 1 && !parent_entry->IsRoot()) {
            // We take the remaining entry from the current node and inserts it into the
            // parent_entry where it replaces (and implicitly deletes) the current node.
            parent_entry->ReplaceNodeWithDataFromEntry(std::move(entries_.begin()->second));
            // WARNING: (this) is deleted here, do not refer to it beyond this point.
        }
    }

    auto& Entries() {
        return entries_;
    }

    const auto& Entries() const {
        return entries_;
    }

    void GetStats(
        PhTreeStats& stats, const EntryT& current_entry, bit_width_t current_depth = 0) const {
        size_t num_children = entries_.size();

        ++stats.n_nodes_;
        ++stats.node_depth_hist_[current_depth];
        ++stats.node_size_log_hist_[32 - detail::CountLeadingZeros(std::uint32_t(num_children))];
        stats.n_total_children_ += num_children;
        stats.q_total_depth_ += current_depth;

        for (auto& entry : entries_) {
            auto& child = entry.second;
            if (child.IsNode()) {
                auto child_infix_len = child.GetNodeInfixLen(current_entry.GetNodePostfixLen());
                ++stats.infix_hist_[child_infix_len];
                auto& sub = child.GetNode();
                sub.GetStats(stats, child, current_depth + 1 + child_infix_len);
            } else {
                ++stats.q_n_post_fix_n_[current_depth];
                ++stats.size_;
            }
        }
    }

    size_t CheckConsistency(const EntryT& current_entry, bit_width_t current_depth = 0) const {
        // Except for a root node if the tree has <2 entries.
        assert(entries_.size() >= 2 || current_depth == 0);
        size_t num_entries_local = 0;
        size_t num_entries_children = 0;
        bool prev_is_node = false;
        hc_pos_t prev_hc_pos = std::numeric_limits<hc_pos_t>::max();
        for (auto& entry : entries_) {
            auto& child = entry.second;
            if (child.IsNode()) {
                assert(entry.first != prev_hc_pos);
                prev_hc_pos = entry.first;
                prev_is_node = true;
                auto& sub = child.GetNode();
                // Check node consistency
                auto sub_infix_len = child.GetNodeInfixLen(current_entry.GetNodePostfixLen());
                assert(
                    sub_infix_len + 1 + child.GetNodePostfixLen() ==
                    current_entry.GetNodePostfixLen());
                num_entries_children +=
                    sub.CheckConsistency(child, current_depth + 1 + sub_infix_len);
            } else {
                assert(prev_is_node == false || entry.first != prev_hc_pos);
                prev_hc_pos = entry.first;
                prev_is_node = false;
                ++num_entries_local;
            }
        }

        // Check node center
        auto post_len = current_entry.GetNodePostfixLen();
        if (post_len == detail::MAX_BIT_WIDTH<SCALAR> - 1) {
            for (auto d : current_entry.GetKey()) {
                assert(d == 0);
            }
        } else {
            for (auto d : current_entry.GetKey()) {
                assert(((d >> post_len) & 0x1) == 1 && "Last bit of node center must be `1`");
                using us_t = std::make_unsigned_t<decltype(d)>;
                us_t d2 = static_cast<us_t>(d);
                assert(((d2 >> post_len) << post_len) == d2 && "postlen bits must all be `0`");
            }
        }

        return num_entries_local + num_entries_children;
    }

  private:
    template <typename... Args>
    auto& WriteValue(hc_pos_t hc_pos, const KeyT& new_key, Args&&... args) {
        return entries_.emplace(hc_pos, EntryT{new_key, std::forward<Args>(args)...})->second;
    }

    void WriteEntry(hc_pos_t hc_pos, EntryT& entry) {
        if (entry.IsNode()) {
            auto postfix_len = entry.GetNodePostfixLen();
            entries_.emplace(hc_pos, EntryT{entry.GetKey(), entry.ExtractNode(), postfix_len});
        } else {
            entries_.emplace(hc_pos, EntryT{entry.GetKey(), entry.ExtractValue()});
        }
    }

    /*
     * Handles the case where we want to insert a new entry into a node but the node already
     * has an entry in that position.
     * @param existing_entry The current entry in the node
     * @param is_inserted Output: This will be set to 'true' by this function if a new entry was
     * inserted by this function.
     * @param new_key The key of the entry to be inserted
     * @param args The constructor arguments for a new value T of a the new entry to be inserted
     * @return A Entry that may contain a child node, a newly created entry or an existing entry.
     * A child node indicates that no entry was inserted, but the caller should try inserting into
     * the child node. A newly created entry (indicated by is_inserted=true) indicates successful
     * insertion. An existing entry (indicated by is_inserted=false) indicates that there is already
     * an entry with the exact same key as new_key, so insertion has failed.
     */
    template <typename... Args>
    auto& HandleCollision(
        EntryT& entry,
        bool& is_inserted,
        const KeyT& new_key,
        bit_width_t current_postfix_len,
        Args&&... args) {
        // We have two entries in the same location (local pos).
        // Now we need to compare the keys, respectively the prefix of the subnode.
        // If they match, we return the entry for further traversal.
        bool is_node = entry.IsNode();
        if (is_node && !entry.HasNodeInfix(current_postfix_len)) {
            // No infix conflict (because infix has length=0), just traverse subnode
            return entry;
        }

        bit_width_t max_conflicting_bits = detail::NumberOfDivergingBits(new_key, entry.GetKey());
        if (is_node && max_conflicting_bits <= entry.GetNodePostfixLen() + 1) {
            // TODO merge with condition above?
            // TODO why is this needed for MMM?
            return entry;
        }
        is_inserted = true;
        return InsertSplit(entry, new_key, max_conflicting_bits, std::forward<Args>(args)...);
    }

    template <typename... Args>
    auto& InsertSplit(
        EntryT& current_entry,
        const KeyT& new_key,
        bit_width_t max_conflicting_bits,
        Args&&... args) {
        // TODO this is not very elegant...  -> handle identical keys
        bit_width_t new_postfix_len = max_conflicting_bits == 0 ? 0 : max_conflicting_bits - 1;
        hc_pos_t pos_sub_1 = detail::CalcPosInArray(new_key, new_postfix_len);
        hc_pos_t pos_sub_2 = detail::CalcPosInArray(current_entry.GetKey(), new_postfix_len);

        // Move key/value into subnode
        Node new_sub_node{};
        new_sub_node.WriteEntry(pos_sub_2, current_entry);
        auto& new_entry = new_sub_node.WriteValue(pos_sub_1, new_key, std::forward<Args>(args)...);

        // Insert new node into local node
        current_entry.SetNode(std::move(new_sub_node), new_postfix_len);
        return new_entry;
    }

    template <typename... Args>
    EntryT& InsertOnSplit(
        EntryT& entry,
        bool& is_inserted,
        const KeyT& new_key,
        bit_width_t postfix_len,
        Args&&... args) {
        // TODO multiple forwarding....?!?!?!?
        EntryT* current =
            &HandleCollision(entry, is_inserted, new_key, postfix_len, std::forward<Args>(args)...);
        while (current->IsNode()) {
            current = &current->GetNode().EmplaceUnchecked(
                is_inserted, new_key, current->GetNodePostfixLen(), std::forward<Args>(args)...);
        }
        assert(is_inserted);
        return *current;
    }

    template <typename... Args>
    auto& SplitQuadrantAndInsert(
        const EntryIterator<DIM, EntryT>& iter_start,
        bool& is_inserted,
        bit_width_t postfix_len,
        const KeyT& new_key,
        EntryT* parent_entry,
        Args&&... args) {
        assert(iter_start != entries_.end());
        auto start = iter_start;
        hc_pos_t hc_pos = start->first;

        // We turn the 'start' into a subnode. Then we move all entries between 'second' and 'end'
        // into 'start'.
        auto current = start;
        ++current;  // start with second entry
        while (current != entries_.end() && current->first == hc_pos) {
            assert(current->first == start->first);
            bool is_moved = false;
            InsertOnSplit(
                start->second,
                is_moved,
                current->second.GetKey(),
                postfix_len,
                std::move(current->second.ExtractValue()));
            assert(is_moved);
            ++current;
        }

        // now before erasing the entries (which may invalidate all iterators) we insert the new
        // entry
        EntryT* result = &InsertOnSplit(
            start->second, is_inserted, new_key, postfix_len, std::forward<Args>(args)...);

        // erase all entries that have been moved
        ++start;
        // TODO check later. THis happens when we have only one entry at hc_pos and turn it
        //   into a sub-node before adding a duplicate (or any key with same hc_pos)
        if (start != current) {
            entries_.erase(start, current);
            CheckUnderflow(parent_entry);
        }

        return *result;
    }

    /*
     * Checks whether an entry's key matches another key. For Non-node entries this simply means
     * comparing the two keys. For entries that contain nodes, we only compare the prefix.
     * @param entry An entry
     * @param key A key
     * @return 'true' iff the relevant part of the key matches (prefix for nodes, whole key for
     * other entries).
     */
    // TODO is this still used for VALUEs?
    bool DoesEntryMatch(
        const EntryT& entry, const KeyT& key, const bit_width_t parent_postfix_len) const {
        if (entry.IsNode()) {
            if (entry.HasNodeInfix(parent_postfix_len)) {
                return detail::KeyEquals(entry.GetKey(), key, entry.GetNodePostfixLen() + 1);
            }
            return true;
        }
        return entry.GetKey() == key;
    }

    EntryMap<DIM, EntryT> entries_;
};

}  // namespace improbable::phtree::v20
#endif  // PHTREE_V20_NODE_H
