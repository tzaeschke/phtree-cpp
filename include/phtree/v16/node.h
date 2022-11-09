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

#ifndef PHTREE_V16_NODE_H
#define PHTREE_V16_NODE_H

#include "entry.h"
#include "phtree/common/b_plus_tree_hash_map.h"
#include "phtree/common/common.h"
#include "phtree_v16.h"
#include <iostream>
#include <map>
#include <set>

namespace improbable::phtree::v16 {

//#define FLEX 1
#define FLEX_SET 1

/*
 * We provide different implementations of the node's internal entry set.
 * All implementations are equivalent to "std::map<hc_pos_t, Entry>" which can be used as
 * a plugin example for verification.
 *
 * - `array_map` is the fastest, but has O(2^DIM) space complexity. This can be very wasteful
 *   because many nodes may have only 2 entries.
 *   Also, iteration depends on some bit operations and is also O(DIM) per step if the CPU/compiler
 *   does not support CTZ (count trailing bits).
 * - `sparse_map` is slower, but requires only O(n) memory (n = number of entries/children).
 *   However, insertion/deletion is O(n), i.e. O(2^DIM) time complexity in the worst case.
 * - 'b_plus_tree_map` is the least efficient for small node sizes but scales best with larger
 *   nodes and dimensionality. Remember that n_max = 2^DIM.
 */

template <typename EntryT>
struct SetEntry {
    hc_pos_t first;
    EntryT second;

    SetEntry() : first{}, second{} {};

    SetEntry(hc_pos_t hc_pos) : first{hc_pos}, second{{}} {};

//    template <typename EntryT2 = EntryT>
//    SetEntry(hc_pos_t hc_pos, EntryT2&& e) : first{hc_pos}, second{std::forward<EntryT>(e)} {};

    template <typename... Args, typename EntryT2 = EntryT>
    explicit SetEntry(hc_pos_t hc_pos, Args&&... args) noexcept
            : first{hc_pos}, second(std::forward<Args>(args)...) {}


//    bool operator==(const SetEntry<EntryT>& rhs) const noexcept {
//        return first == rhs.first && second.GetKey() == rhs.second.GetKey();
//    }
};
template <typename EntryT>
struct equal_to_content {
    bool operator()(
            const SetEntry<EntryT>& x1,
            const SetEntry<EntryT>& x2) const {
        return (*x1) == (*x2);
    }
};
//template <typename EntryT>
//struct less_content {
//    bool operator()(
//            const SetEntry<EntryT>& x1,
//            const SetEntry<EntryT>& x2) const {
//        return x1.first < x2.first || x1.first == x2.first && x1.second.GetKey()[0] < x2.second.GetKey()[0];
//    }
//};
}  // namespace phtree_multimap_d_test_filter

//namespace std {
//    template <>
//    struct hash<improbable::phtree::v16::SetEntry> {
//        size_t operator()(const phtree_multimap_d_test_filter::Id& x) const {
//            return std::hash<int>{}(x._i);
//        }
//    };
//};  // namespace std

namespace improbable::phtree::v16 {

    template <typename EntryT>
    struct less_content {
        bool operator()(
                const SetEntry<EntryT>& x1,
                const SetEntry<EntryT>& x2) const {
            auto& e1 = x1.second;
            auto& e2 = x2.second;
            if (e1.IsValue() && e2.IsValue()) {
                return x1.first < x2.first || x1.first == x2.first && x1.second.GetKey() < x2.second.GetKey();
            }
            using SCALAR = std::remove_reference_t<decltype(e1.GetKey()[0])>;
            const auto mask1 = MAX_MASK<SCALAR> << (e1.IsNode() ? (e1.GetNodePostfixLen() + 1) : 0);
            const auto mask2 = MAX_MASK<SCALAR> << (e2.IsNode() ? (e2.GetNodePostfixLen() + 1) : 0);
            return KeyLess(e1.GetKey(), e2.GetKey(), mask1 & mask2);
        }
    };

template <dimension_t DIM, typename Entry>
// using EntryMap = typename std::conditional<
//     DIM <= 3,
//     array_map<Entry, (hc_pos_t(1) << DIM)>,
//     typename std::
//         conditional<DIM <= 8, sparse_map<Entry>, b_plus_tree_map<Entry, (hc_pos_t(1) << DIM)>>::
//             type>::type;
// using EntryMap = typename std::conditional<
//    DIM <= 8,
//    b_plus_tree_hash_map<hc_pos_t, Entry>,
//    std::multimap<hc_pos_t, Entry>>::type;

#ifdef FLEX
using EntryMap = std::multimap<hc_pos_t, Entry>;
#elif defined(FLEX_SET)
using EntryMap = std::set<SetEntry<Entry>, less_content<Entry>>;
#else
 using EntryMap = typename std::conditional<
     DIM <= 3,
     array_map<Entry, (hc_pos_t(1) << DIM)>,
     typename std::
         conditional<DIM <= 8, sparse_map<Entry>, b_plus_tree_map<Entry, (hc_pos_t(1) << DIM)>>::
             type>::type;
#endif

template <dimension_t DIM, typename Entry>
using EntryIterator = decltype(EntryMap<DIM, Entry>().begin());
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
#ifdef FLEX
    static constexpr hc_pos_t MAX_SIZE = std::max(hc_pos_t(8), hc_pos_t(1) << DIM);
#endif
    using KeyT = PhPoint<DIM, SCALAR>;
    using EntryT = Entry<DIM, T, SCALAR>;

  public:
    Node() : entries_{} {}

    // Nodes should never be copied!
    Node(const Node&) = delete;
    Node(Node&&) = delete;
    Node& operator=(const Node&) = delete;
    Node& operator=(Node&&) = delete;

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
    EntryT& Emplace(bool& is_inserted, const KeyT& key, bit_width_t postfix_len, Args&&... args) {
        hc_pos_t hc_pos = CalcPosInArray(key, postfix_len);
#ifdef FLEX
        size_t match_count = 0;
        {
            // subnode? return for further traversal
            auto iter = entries_.lower_bound(hc_pos);
            if (iter != entries_.end() && iter->second.IsNode() && iter->first == hc_pos) {
                return HandleCollision(
                    iter->second, is_inserted, key, postfix_len, std::forward<Args>(args)...);
            }

            // Does entry exist? -> return as failed insertion
            // TODO Implement as native multimap??? -> What about MAX_SIZE?
            while (iter != entries_.end() && iter->first == hc_pos) {
                ++match_count;
                if (iter->second.GetKey() == key) {
                    // Entry exists ...
                    return iter->second;
                }
                ++iter;
            }
        }

        // split node if it is too large
        // TODO std::cerr << "hhhhhh-------- " << entries_.size() << " >= " << MAX_SIZE << "  DIM=" << DIM << std::endl;
        // TODO DIM < 32!!!
        //if (DIM <32 && entries_.size() >= MAX_SIZE) {
        if (match_count >= 4) {
            Split(postfix_len);
        }

        // search again
        // TODO avoid this, we can e.g. calculate the new iter position.
        auto iter = entries_.lower_bound(hc_pos);
        // Subnode? return for further traversal
        if (iter != entries_.end() && iter->second.IsNode() && iter->first == hc_pos) {
            return HandleCollision(
                iter->second, is_inserted, key, postfix_len, std::forward<Args>(args)...);
        }

        // Does not exist
        is_inserted = true;
        auto entry_iter =
            entries_.emplace_hint(iter, hc_pos, EntryT{key, std::forward<Args>(args)...});
        return entry_iter->second;
#elif defined(FLEX_SET)
        size_t match_count = 0;
        {
            // subnode? return for further traversal
            auto iter = entries_.lower_bound(wrap(hc_pos));
            if (iter != entries_.end() && iter->second.IsNode() && iter->first == hc_pos) {
                EntryT& e = const_cast<EntryT&>(iter->second);  // TODO why cast???
                return HandleCollision(
                        e, is_inserted, key, postfix_len, std::forward<Args>(args)...);
            }

            // Does entry exist? -> return as failed insertion
            // TODO Implement as native multimap??? -> What about MAX_SIZE?
            while (iter != entries_.end() && iter->first == hc_pos) {
                ++match_count;
                if (iter->second.GetKey() == key) {
                    // Entry exists ...
                    return const_cast<EntryT&>(iter->second);  // TODO why cast???
                }
                ++iter;
            }
        }

        // split node if it is too large
        // TODO std::cerr << "hhhhhh-------- " << entries_.size() << " >= " << MAX_SIZE << "  DIM=" << DIM << std::endl;
        // TODO DIM < 32!!!
        //if (DIM <32 && entries_.size() >= MAX_SIZE) {
        if (match_count >= 4) {
            Split(postfix_len);
        }

        // search again
        // TODO avoid this, we can e.g. calculate the new iter position. OR WE WE DIDN'T SPLIT!
        auto iter = entries_.lower_bound(wrap(hc_pos));
        // Subnode? return for further traversal
        if (iter != entries_.end() && iter->second.IsNode() && iter->first == hc_pos) {
            EntryT& e = const_cast<EntryT&>(iter->second);  // TODO why cast???
            return HandleCollision(
                    e, is_inserted, key, postfix_len, std::forward<Args>(args)...);
        }

        // Does not exist
//        is_inserted = true;
//        auto entry_iter =
//                entries_.emplace_hint(iter, hc_pos, EntryT{key, std::forward<Args>(args)...});
//        return const_cast<EntryT&>(entry_iter->second);  // TODO why cast???
        auto result = entries_.emplace(hc_pos, EntryT{key, std::forward<Args>(args)...});
        if (result.second) {
            is_inserted = true;
        }
        return const_cast<EntryT&>(result.first->second);
#else
        auto emplace_result = entries_.try_emplace(hc_pos, key, std::forward<Args>(args)...);
        auto& entry = emplace_result.first->second;
        // Return if emplace succeed, i.e. there was no entry.
        if (emplace_result.second) {
            is_inserted = true;
            return entry;
        }
        return HandleCollision(entry, is_inserted, key, postfix_len, std::forward<Args>(args)...);
#endif
    }

        static SetEntry<EntryT> wrap(hc_pos_t hc_pos) {
            return SetEntry<EntryT>{hc_pos};
        }

        static SetEntry<EntryT> wrap(hc_pos_t hc_pos, const KeyT& key) {
            return SetEntry<EntryT>{hc_pos, key};
        }

        auto lower_bound(hc_pos_t hc_pos) {
            return entries_.lower_bound(wrap(hc_pos));
        }

        auto lower_bound(hc_pos_t hc_pos) const {
            return entries_.lower_bound(wrap(hc_pos));
        }

        /*
     * Returns the value (T or Node) if the entry exists and matches the key. Child nodes are
     * _not_ traversed.
     * @param key The key of the entry
     * @param parent The parent node
     * @return The sub node or null.
     */
    const EntryT* Find(const KeyT& key, bit_width_t postfix_len) const {
        hc_pos_t hc_pos = CalcPosInArray(key, postfix_len);
#if defined(FLEX) || defined(FLEX_SET)
        auto iter = entries_.lower_bound(wrap(hc_pos, key));
        while (iter != entries_.end() && iter->first == hc_pos) {
            if (DoesEntryMatch(iter->second, key, postfix_len)) {
                return &iter->second;
            }
            ++iter;
        }
        return nullptr;
#else
        const auto iter = entries_.find(hc_pos);
        if (iter != entries_.end() && DoesEntryMatch(iter->second, key, postfix_len)) {
            return &iter->second;
        }
        return nullptr;
#endif
    }

    EntryT* Find(const KeyT& key, bit_width_t postfix_len) {
        return const_cast<EntryT*>(static_cast<const Node*>(this)->Find(key, postfix_len));
    }

    EntryIteratorC<DIM, EntryT> FindPrefix(
        const KeyT& prefix, bit_width_t prefix_post_len, bit_width_t node_postfix_len) const {
        assert(prefix_post_len <= node_postfix_len);
        hc_pos_t hc_pos = CalcPosInArray(prefix, node_postfix_len);
        const auto iter = entries_.find(wrap(hc_pos));
        if (iter == entries_.end() || iter->second.IsValue() ||
            iter->second.GetNodePostfixLen() < prefix_post_len) {
            // We compare the infix only if it lies fully within the prefix.
            return entries_.end();
        }

        if (DoesEntryMatch(iter->second, prefix, node_postfix_len)) {
            return {iter};
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
    EntryT* Erase(const KeyT& key, EntryT* parent_entry, bool allow_move_into_parent, bool& found) {
        auto postfix_len = parent_entry->GetNodePostfixLen();
        hc_pos_t hc_pos = CalcPosInArray(key, postfix_len);
#if defined(FLEX) || defined(FLEX_SET)
        auto iter = entries_.lower_bound(wrap(hc_pos));
        while (iter != entries_.end() && iter->first == hc_pos) {
            if (DoesEntryMatch(iter->second, key, postfix_len)) {
                if (iter->second.IsNode()) {
                    return & const_cast<EntryT&>(iter->second);
                }
                entries_.erase(iter);

                found = true;
                if (allow_move_into_parent && GetEntryCount() == 1) {
                    // We take the remaining entry from the current node and inserts it into the
                    // parent_entry where it replaces (and implicitly deletes) the current node.
                    parent_entry->ReplaceNodeWithDataFromEntry(std::move(const_cast<EntryT&>(entries_.begin()->second)));
                    // WARNING: (this) is deleted here, do not refer to it beyond this point.
                }
                return nullptr;
            }
            ++iter;
        }
        return nullptr;
#else
        auto it = entries_.find(hc_pos);
        if (it != entries_.end() && DoesEntryMatch(it->second, key, postfix_len)) {
            if (it->second.IsNode()) {
                return &it->second;
            }
            entries_.erase(it);

            found = true;
            if (allow_move_into_parent && GetEntryCount() == 1) {
                // We take the remaining entry from the current node and inserts it into the
                // parent_entry where it replaces (and implicitly deletes) the current node.
                parent_entry->ReplaceNodeWithDataFromEntry(std::move(entries_.begin()->second));
                // WARNING: (this) is deleted here, do not refer to it beyond this point.
            }
        }
        return nullptr;
#endif
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
        ++stats.node_size_log_hist_[32 - CountLeadingZeros(std::uint32_t(num_children))];
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
        return num_entries_local + num_entries_children;
    }

  private:
    template <typename... Args>
    EntryT& WriteValue(hc_pos_t hc_pos, const KeyT& new_key, Args&&... args) {
#if defined(FLEX)
        return entries_.emplace(hc_pos, EntryT{new_key, std::forward<Args>(args)...})->second;
#elif defined(FLEX_SET)
        return const_cast<EntryT&>(entries_.emplace(hc_pos, new_key, std::forward<Args>(args)...).first->second);
#else
        return entries_.try_emplace(hc_pos, new_key, std::forward<Args>(args)...).first->second;
#endif
    }

    void WriteEntry(hc_pos_t hc_pos, EntryT& entry) {
#if defined(FLEX) || defined(FLEX_SET)
        if (entry.IsNode()) {
            auto postfix_len = entry.GetNodePostfixLen();
            entries_.emplace(hc_pos, EntryT{entry.GetKey(), entry.ExtractNode(), postfix_len});
        } else {
            entries_.emplace(hc_pos, EntryT{entry.GetKey(), entry.ExtractValue()});
        }
#else
        if (entry.IsNode()) {
            auto postfix_len = entry.GetNodePostfixLen();
            entries_.try_emplace(hc_pos, entry.GetKey(), entry.ExtractNode(), postfix_len);
        } else {
            entries_.try_emplace(hc_pos, entry.GetKey(), entry.ExtractValue());
        }
#endif
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
    EntryT& HandleCollision(
        EntryT& existing_entry,
        bool& is_inserted,
        const KeyT& new_key,
        bit_width_t current_postfix_len,
        Args&&... args) {
        assert(!is_inserted);
        // We have two entries in the same location (local pos).
        // Now we need to compare the keys.
        // If they are identical, we simply return the entry for further traversal.
        if (existing_entry.IsNode()) {
            if (existing_entry.HasNodeInfix(current_postfix_len)) {
                bit_width_t max_conflicting_bits =
                    NumberOfDivergingBits(new_key, existing_entry.GetKey());
                if (max_conflicting_bits > existing_entry.GetNodePostfixLen() + 1) {
                    is_inserted = true;
                    return InsertSplit(
                        existing_entry, new_key, max_conflicting_bits, std::forward<Args>(args)...);
                }
            }
            // No infix conflict, just traverse subnode
        } else {
            bit_width_t max_conflicting_bits =
                NumberOfDivergingBits(new_key, existing_entry.GetKey());
            if (max_conflicting_bits > 0) {
                is_inserted = true;
                return InsertSplit(
                    existing_entry, new_key, max_conflicting_bits, std::forward<Args>(args)...);
            }
            // perfect match -> return existing
        }
        return existing_entry;
    }

    template <typename... Args>
    EntryT& InsertSplit(
        EntryT& current_entry,
        const KeyT& new_key,
        bit_width_t max_conflicting_bits,
        Args&&... args) {
        bit_width_t new_postfix_len = max_conflicting_bits - 1;
        auto new_sub_node = std::make_unique<Node>();
        hc_pos_t pos_sub_1 = CalcPosInArray(new_key, new_postfix_len);
        hc_pos_t pos_sub_2 = CalcPosInArray(current_entry.GetKey(), new_postfix_len);

        // Move key/value into subnode
        new_sub_node->WriteEntry(pos_sub_2, current_entry);
        auto& new_entry = new_sub_node->WriteValue(pos_sub_1, new_key, std::forward<Args>(args)...);

        // Insert new node into local node
        current_entry.SetNode(std::move(new_sub_node), new_postfix_len);
        return new_entry;
    }

#if defined(FLEX) || defined(FLEX_SET)
    void Split(bit_width_t postfix_len) {
        auto iter = entries_.begin();
        int nIter = 0;
        hc_pos_t current_hc_pos = std::numeric_limits<hc_pos_t>::max();

        // find start
        auto start = iter;
        int nStart = nIter;
        while (iter != entries_.end() && iter->first != current_hc_pos) {
            current_hc_pos = iter->first;
            start = iter;
            nStart = nIter;
            ++iter;
            ++nIter;
        }

        // find "end"
        auto end = iter;
        int nEnd = nIter;
        while (iter != entries_.end() && iter->first == current_hc_pos) {
            current_hc_pos = iter->first;
            ++iter;
            ++nIter;
            end = iter;
            nEnd = nIter;
        }

        // insert "second" entry into "first"
        auto second = start;
        ++second;
        int nSecond = nStart;
        ++nSecond;
        bool dummy;

        auto current = second;
        int nCurrent = nSecond;
//        std::cout << "iter= " << nIter << " current= " << nCurrent << " second=" << nSecond
//                  << " start= " << nStart << " end= " << nEnd << " size=" << entries_.size() << "  hcpos=" << start->first << std::endl;
        while (current != end) {
            assert(current->first == start->first);
            bool is_inserted = false;
 //           std::cout << " current= " << nCurrent << "  key=" << current->second.GetKey() << std::endl;
            HandleCollision(
                const_cast<EntryT&>(start->second),                            // EntryT& existing_entry,
                is_inserted,                                    // bool& is_inserted,
                current->second.GetKey(),                  // const KeyT& new_key,
                postfix_len,                              // bit_width_t current_postfix_len,
                const_cast<EntryT&>(current->second).ExtractValue()  // Args&&... args
            );
            // If it was not inserted via collision handling (=collision with sub-node prefix)
            // then we insert it directly into the new subnode.
            if (!is_inserted) {
                auto* current_entry = &start->second;
                while (current_entry->IsNode()) {
//                    std::cout << "  current= " << nCurrent << std::endl;
                    current_entry = &current_entry->GetNode().Emplace(
                        is_inserted,
                        current->second.GetKey(),
                        current_entry->GetNodePostfixLen(),
                        const_cast<EntryT&>(current->second).ExtractValue());
                }
            }
            // TODO we could have a nested overflow....
            // assert(is_inserted);
            ++current;
            ++nCurrent;
        }

//        std::cout << "Before" << std::endl;
//        for (auto& x : entries_) {
//            std::cout << "   e=" << x.first << " -> " << x.second.GetKey() << std::endl;
//        }
        entries_.erase(second, end);
//        std::cout << "After" << std::endl;
//        for (auto& x : entries_) {
//            std::cout << "   e=" << x.first << " -> " << x.second.GetKey() << std::endl;
//        }
        auto next = start;
        ++next;
        assert(next == entries_.end() || start->first < next->first);
        //assert((start++) == entries_.end() || start->first < (++start)->first);
    }
#endif

    /*
     * Checks whether an entry's key matches another key. For Non-node entries this simply means
     * comparing the two keys. For entries that contain nodes, we only compare the prefix.
     * @param entry An entry
     * @param key A key
     * @return 'true' iff the relevant part of the key matches (prefix for nodes, whole key for
     * other entries).
     */
    bool DoesEntryMatch(
        const EntryT& entry, const KeyT& key, const bit_width_t parent_postfix_len) const {
        if (entry.IsNode()) {
            if (entry.HasNodeInfix(parent_postfix_len)) {
                const bit_mask_t<SCALAR> mask = MAX_MASK<SCALAR> << (entry.GetNodePostfixLen() + 1);
                return KeyEquals(entry.GetKey(), key, mask);
            }
            return true;
        }
        return entry.GetKey() == key;
    }

    EntryMap<DIM, EntryT> entries_;
};

}  // namespace improbable::phtree::v16
#endif  // PHTREE_V16_NODE_H
