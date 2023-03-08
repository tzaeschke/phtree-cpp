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

#ifndef PHTREE_V20_PHTREE_V20_H
#define PHTREE_V20_PHTREE_V20_H

#include "debug_helper_v20.h"
#include "for_each.h"
#include "for_each_hc.h"
#include "iterator_full.h"
#include "iterator_hc.h"
#include "iterator_knn_hs.h"
#include "iterator_with_key.h"
#include "iterator_with_parent.h"
#include "node.h"

namespace improbable::phtree::v20 {

/*
 * The PH-Tree is an ordered index on an n-dimensional space (quad-/oct-/2^n-tree) where each
 * dimension is (by default) indexed by a 64 bit integer. The index ordered follows z-order / Morton
 * order. The index is effectively a 'map', i.e. each key is associated with at most one value.
 *
 * Keys are points in n-dimensional space.
 *
 * This API behaves similar to std::map, see function descriptions for details.
 *
 * Loosely based on PH-Tree Java, V20, see http://www.phtree.org
 *
 * See also :
 * - T. Zaeschke, C. Zimmerli, M.C. Norrie:
 *   "The PH-Tree -- A Space-Efficient Storage Structure and Multi-Dimensional Index", (SIGMOD 2014)
 * - T. Zaeschke: "The PH-Tree Revisited", (2015)
 * - T. Zaeschke, M.C. Norrie: "Efficient Z-Ordered Traversal of Hypercube Indexes" (BTW 2017).
 *
 * @tparam T            Value type.
 * @tparam DIM          Dimensionality. This is the number of dimensions of the space to index.
 * @tparam CONVERT      A converter class with a 'pre()' and a 'post()' function. 'pre()' translates
 *                      external KEYs into the internal PhPoint<DIM, SCALAR_BITS> type. 'post()'
 *                      translates the PhPoint<DIM, SCALAR_BITS> back to the external KEY type.
 */
template <dimension_t DIM, typename T, typename CONVERT = ConverterNoOp<DIM, scalar_64_t>>
class PhTreeV20 {
    friend PhTreeDebugHelper;
    using bit_width_t = detail::bit_width_t;
    using ScalarExternal = typename CONVERT::ScalarExternal;
    using ScalarInternal = typename CONVERT::ScalarInternal;
    using KeyT = typename CONVERT::KeyInternal;
    using EntryT = Entry<DIM, T, ScalarInternal>;
    using NodeT = Node<DIM, T, ScalarInternal>;

  public:
    static_assert(!std::is_reference<T>::value, "Reference type value are not supported.");
    static_assert(std::is_signed<ScalarInternal>::value, "ScalarInternal must be a signed type");
    static_assert(
        std::is_integral<ScalarInternal>::value, "ScalarInternal must be an integral type");
    static_assert(
        std::is_arithmetic<ScalarExternal>::value, "ScalarExternal must be an arithmetic type");
    static_assert(DIM >= 1 && DIM <= 63, "This PH-Tree supports between 1 and 63 dimensions");

    explicit PhTreeV20(CONVERT* converter)
    : num_entries_{0}
    , root_{{}, NodeT{}, detail::MAX_BIT_WIDTH<ScalarInternal> - 1}
    , converter_{converter} {}

    PhTreeV20(const PhTreeV20& other) = delete;
    PhTreeV20& operator=(const PhTreeV20& other) = delete;
    PhTreeV20(PhTreeV20&& other) noexcept = default;
    PhTreeV20& operator=(PhTreeV20&& other) noexcept = default;
    ~PhTreeV20() noexcept = default;

    /*
     *  Attempts to build and insert a key and a value into the tree.
     *
     *  @param key The key for the new entry.
     *
     *  @param args  Arguments used to generate a new value.
     *
     *  @return  A pair, whose first element points  to the possibly inserted pair,
     *           and whose second element is a bool that is true if the pair was actually inserted.
     *
     * This function attempts to build and insert a (key, value) pair into the tree. The PH-Tree is
     * effectively a map, so if an entry with the same key was already in the tree, returns that
     * entry instead of inserting a new one.
     */
    template <typename... Args>
    std::pair<T&, bool> try_emplace(const KeyT& key, Args&&... args) {
        auto* entry = &root_;
        bool is_inserted = false;
        while (entry->IsNode()) {
            entry = &entry->GetNode().Emplace(is_inserted, key, entry, std::forward<Args>(args)...);
        }
        num_entries_ += is_inserted;
        return {entry->GetValue(), is_inserted};
    }

    /*
     * The try_emplace(hint, key, value) method uses an iterator as hint for insertion.
     * The hint is ignored if it is not useful or is equal to end().
     *
     * Iterators should normally not be used after the tree has been modified. As an exception to
     * this rule, an iterator can be used as hint if it was previously used with at most one call
     * to erase() and if no other modifications occurred.
     * The following is valid:
     *
     * // Move value from key1 to key2
     * auto iter = tree.find(key1);
     * auto value = iter.second(); // The value may become invalid in erase()
     * erase(iter);
     * try_emplace(iter, key2, value);  // the iterator can still be used as hint here
     */
    template <typename ITERATOR, typename... Args>
    std::pair<T&, bool> try_emplace(const ITERATOR& iterator, const KeyT& key, Args&&... args) {
        if constexpr (!std::is_same_v<ITERATOR, IteratorWithParent<T, CONVERT>>) {
            return try_emplace(key, std::forward<Args>(args)...);
        } else {
            // This function can be used to insert a value close to a known value
            // or close to a recently removed value. The hint can only be used if the new key is
            // inside one of the nodes provided by the hint iterator.
            // The idea behind using the 'parent' is twofold:
            // - The 'parent' node is one level above the iterator position, it is spatially
            //   larger and has a better probability of containing the new position, allowing for
            //   fast track try_emplace.
            // - Using 'parent' allows a scenario where the iterator was previously used with
            //   erase(iterator). This is safe because erase() will never erase the 'parent' node.

            if (!iterator.__GetParentNodeEntry()) {
                // No hint available, use standard try_emplace()
                return try_emplace(key, std::forward<Args>(args)...);
            }

            auto* parent_entry = iterator.__GetParentNodeEntry();
            if (detail::NumberOfDivergingBits(key, parent_entry->GetKey()) >
                parent_entry->GetNodePostfixLen() + 1) {
                // replace higher up in the tree
                return try_emplace(key, std::forward<Args>(args)...);
            }

            // replace in node
            auto* entry = parent_entry;
            bool is_inserted = false;
            while (entry->IsNode()) {
                entry =
                    &entry->GetNode().Emplace(is_inserted, key, entry, std::forward<Args>(args)...);
            }
            num_entries_ += is_inserted;
            return {entry->GetValue(), is_inserted};
        }
    }

    /*
     * See std::map::insert().
     *
     * @return a pair consisting of the inserted element (or to the element that prevented the
     * insertion) and a bool denoting whether the insertion took place.
     */
    std::pair<T&, bool> insert(const KeyT& key, const T& value) {
        return try_emplace(key, value);
    }

    /*
     * @return the value stored at position 'key'. If no such value exists, one is added to the tree
     * and returned.
     */
    T& operator[](const KeyT& key) {
        return try_emplace(key).first;
    }

    /*
     * Analogous to map:count().
     *
     * @return '1', if a value is associated with the provided key, otherwise '0'.
     */
    size_t count(const KeyT& key) const {
        if (empty()) {
            return 0;
        }
        auto* current_node = &root_;
        //        while (current_entry && current_entry->IsNode()) {
        //            current_entry = current_entry->GetNode().FindC(key,
        //            current_entry->GetNodePostfixLen());
        //        }
        bool is_found = true;  // Yep! Otherwise we could use a do-while loop. // TODO ?
        EntryIteratorC<DIM, EntryT> iter =
            root_.GetNode().LowerBound(key, root_.GetNodePostfixLen(), is_found);
        ;
        // TODO use iter.IsEnd() i.o. is_found ?!?!?
        while (is_found && iter->second.IsNode()) {
            current_node = &iter->second;
            iter =
                iter->second.GetNode().LowerBound(key, iter->second.GetNodePostfixLen(), is_found);
        }
        int n = 0;
        // TODO check for hc_pos change to abort early
        while (iter != current_node->GetNode().End()) {
            if (iter->second.GetKey() == key) {
                ++n;
            }
            ++iter;
        }
        return n;
    }

    /*
     * Analogous to map:find().
     *
     * Get an entry associated with a k dimensional key.
     * @param key the key to look up
     * @return an iterator that points either to the associated value or to {@code end()} if the key
     * was found
     */
    auto find(const KeyT& key) const {
        auto* current_node = &root_;
        bool is_found = true;
        EntryIteratorC<DIM, EntryT> iter =
            root_.GetNode().LowerBound(key, root_.GetNodePostfixLen(), is_found);
        ;
        while (is_found && iter->second.IsNode()) {
            current_node = &iter->second;
            iter =
                iter->second.GetNode().LowerBound(key, iter->second.GetNodePostfixLen(), is_found);
        }
        if (is_found) {  // iter != current_node->GetNode().End()
            return IteratorWithKey<T, CONVERT>(key, iter, &iter->second, current_node, converter_);
        }
        return IteratorWithKey<T, CONVERT>(key, {}, nullptr, nullptr, converter_);
    }

    auto find(const KeyT& key, const T& value) const {
        auto* current_node = &root_;
        bool is_found = true;
        EntryIteratorC<DIM, EntryT> iter =
            root_.GetNode().LowerBound(key, value, root_.GetNodePostfixLen(), is_found);
        while (is_found && iter->second.IsNode()) {
            current_node = &iter->second;
            iter = iter->second.GetNode().LowerBound(
                key, value, iter->second.GetNodePostfixLen(), is_found);
        }
        if (is_found) {  // iter != current_node->GetNode().End()
            // TODO parent node?
            return IteratorWithParent<T, CONVERT>(&iter->second, current_node, nullptr, converter_);
        }
        return IteratorWithParent<T, CONVERT>(nullptr, nullptr, nullptr, converter_);
    }

    /*
     * See std::map::erase(). Removes any value associated with the provided key.
     *
     * @return '1' if a value was found, otherwise '0'.
     */
    // TODO is this method still needed (i.e. without "value") ???
    //    size_t erase2(const KeyT& key) {
    //        auto* entry = &root_;
    //        // We do not want the root entry to be modified. The reason is simply that a lot of
    //        the
    //        // code in this class becomes simpler if we can assume the root entry to contain a
    //        node. bool found = false; while (entry) {
    //            entry = entry->GetNode().Erase(key, entry, found);
    //        }
    //        num_entries_ -= found;
    //        return found;
    //    }

    size_t erase(const KeyT& key, const T& value) {
        auto* entry = &root_;
        // We do not want the root entry to be modified. The reason is simply that a lot of the
        // code in this class becomes simpler if we can assume the root entry to contain a node.
        size_t found = 0;
        while (entry) {
            entry = entry->GetNode().Erase(key, value, entry, found);
        }
        num_entries_ -= found;
        return found;
    }

    /*
     * See std::map::erase(). Removes any value at the given iterator location.
     *
     * WARNING
     * While this is guaranteed to work correctly, only iterators returned from find()
     * will result in erase(iterator) being faster than erase(key).
     * Iterators returned from other functions may be optimized in a future version.
     *
     * @return '1' if a value was found, otherwise '0'.
     */
    template <typename ITERATOR>
    size_t erase(const ITERATOR& iterator) {
        if (iterator.IsEnd()) {
            return 0;
        }
        // TODO this needs cleaning up...?!?!? -> Add node-iterator to IteratorWithParent
        //   We should return an iterable iterator here. IterWIthParent is _not_ iterable!
        if constexpr (std::is_same_v<ITERATOR, IteratorWithParent<T, CONVERT>>) {
            const auto& iter_rich = static_cast<const IteratorWithParent<T, CONVERT>&>(iterator);
            if (!iter_rich.__GetNodeEntry() || iter_rich.__GetNodeEntry() == &root_) {
                // Do _not_ use the root entry, see erase(key). Start searching from the top.
                return erase(iter_rich.GetEntry()->GetKey(), iterator.GetEntry()->GetValue());
            }
            size_t found = 0;
            EntryT* entry = iter_rich.__GetNodeEntry();
            // The loop is a safeguard for find_two_mm which may return slightly wrong iterators.
            while (entry != nullptr) {
                entry = entry->GetNode().Erase(
                    iter_rich.GetEntry()->GetKey(), iter_rich.GetEntry()->GetValue(), entry, found);
            }
            num_entries_ -= found;
            return found;
        }
        // There may be no entry because not every iterator sets it.
        return erase(iterator.GetEntry()->GetKey(), iterator.GetEntry()->GetValue());
    }

    /*
     * Relocate (move) ONE entry from one position to another, subject to a predicate.
     *
     * @param old_key
     * @param new_key
     * @param predicate
     *
     * @return  A pair, whose first element points to the possibly relocated value, and
     *          whose second element is a bool that is true if the value was actually relocated.
     */
    template <typename PRED>
    size_t relocate_if(const KeyT& old_key, const KeyT& new_key, PRED&& pred) {
        bit_width_t n_diverging_bits = detail::NumberOfDivergingBits(old_key, new_key);

        EntryT* current_entry = &root_;           // An entry.
        EntryT* old_node_entry = nullptr;         // Node that contains entry to be removed
        EntryT* new_node_entry_parent = nullptr;  // Parent of the old_node_entry
        EntryT* new_node_entry = nullptr;         // Node that will contain  new entry

        auto iter = root_.GetNode().FindIter(old_key, root_.GetNodePostfixLen());
        // Find node for removal
        while (current_entry && current_entry->IsNode()) {
            old_node_entry = current_entry;
            auto postfix_len = old_node_entry->GetNodePostfixLen();
            if (postfix_len + 1 >= n_diverging_bits) {
                new_node_entry_parent = new_node_entry;
                new_node_entry = old_node_entry;
            }
            // TODO stop earlier, we are going to have to redo this after insert....
            iter = current_entry->GetNode().FindIter(old_key, postfix_len);
            current_entry = iter == current_entry->GetNode().End() ? nullptr : &iter->second;
        }
        auto iter_end = old_node_entry->GetNode().End();
        EntryT* old_entry = current_entry;  // Entry to be removed

        // TODO abort early for hc_pos mismatch.
        //  Maybe we can return an end-iterator from Find()? E.g. lowerbound hc_pos + 1?
        while (iter != iter_end) {
            old_entry = &iter->second;  // TODO use ref& ?
            if (old_entry->GetKey() == old_key && pred(old_entry->GetValue())) {
                // TODO skip this comparison for first entry? Or don't do it inside FInd()?
                break;
            }
            ++iter;
        }
        if (iter == iter_end) {
            // No matching entry found.
            return 0;
        }

        // Are the keys equal? -> same entry
        if (n_diverging_bits == 0) {
            return 1;
        }
        // Is the quadrant the same? -> same entry
        if (old_node_entry->GetNodePostfixLen() >= n_diverging_bits) {
            old_entry->SetKey(new_key);
            return 1;
        }

        // Extract value
        T value = old_entry->ExtractValue();

        // Erase old entry
        // We do not check for underflow if we are reinserting into the same node
        bool check_underflow = old_node_entry != new_node_entry;
        old_node_entry->GetNode().Erase(iter, old_node_entry, check_underflow);

        // insert
        bool is_inserted = false;
        auto new_entry = new_node_entry;
        while (new_entry->IsNode()) {
            // TODO multiple move()....?!?!?
            new_entry =
                &new_entry->GetNode().Emplace(is_inserted, new_key, new_entry, std::move(value));
        }
        assert(is_inserted);  // TODO remove bool flag?

        if (!check_underflow && old_node_entry->IsNode()) {
            // If we skipped the check above, we have to check now...
            old_node_entry->GetNode().CheckUnderflow(old_node_entry);
        }

        return 1;
    }

    /*
     * Iterates over all entries in the tree. The optional filter allows filtering entries and nodes
     * (=sub-trees) before returning / traversing them. By default, all entries are returned. Filter
     * functions must implement the same signature as the default 'FilterNoOp'.
     *
     * @param callback The callback function to be called for every entry that matches the query.
     * The callback requires the following signature: callback(const PhPointD<DIM> &, const T &)
     * @param filter An optional filter function. The filter function allows filtering entries and
     * sub-nodes before they are returned or traversed. Any filter function must follow the
     * signature of the default 'FilterNoOp`.
     */
    template <typename CALLBACK, typename FILTER = FilterNoOp>
    void for_each(CALLBACK&& callback, FILTER&& filter = FILTER()) {
        ForEach<T, CONVERT, CALLBACK, FILTER>(
            converter_, std::forward<CALLBACK>(callback), std::forward<FILTER>(filter))
            .Traverse(root_);
    }

    template <typename CALLBACK, typename FILTER = FilterNoOp>
    void for_each(CALLBACK&& callback, FILTER&& filter = FILTER()) const {
        ForEach<T, CONVERT, CALLBACK, FILTER>(
            converter_, std::forward<CALLBACK>(callback), std::forward<FILTER>(filter))
            .Traverse(root_);
    }

    /*
     * Performs a rectangular window query. The parameters are the min and max keys which
     * contain the minimum respectively the maximum keys in every dimension.
     * @param query_box The query window.
     * @param callback The callback function to be called for every entry that matches the query.
     * The callback requires the following signature: callback(const PhPoint<DIM> &, const T &)
     * @param filter An optional filter function. The filter function allows filtering entries and
     * sub-nodes before they are returned or traversed. Any filter function must follow the
     * signature of the default 'FilterNoOp`.
     */
    template <typename CALLBACK, typename FILTER = FilterNoOp>
    void for_each(
        // TODO check copy elision
        const PhBox<DIM, ScalarInternal> query_box,
        CALLBACK&& callback,
        FILTER&& filter = FILTER()) const {
        auto pair = find_starting_node(query_box);
        ForEachHC<T, CONVERT, CALLBACK, FILTER>(
            query_box.min(),
            query_box.max(),
            converter_,
            std::forward<CALLBACK>(callback),
            std::forward<FILTER>(filter))
            .Traverse(*pair.first, &pair.second);
    }

    /*
     * Iterates over all entries in the tree. The optional filter allows filtering entries and nodes
     * (=sub-trees) before returning / traversing them. By default all entries are returned. Filter
     * functions must implement the same signature as the default 'FilterNoOp'.
     *
     * @return an iterator over all (filtered) entries in the tree,
     */
    template <typename FILTER = FilterNoOp>
    auto begin(FILTER&& filter = FILTER()) const {
        return IteratorFull<T, CONVERT, FILTER>(root_, converter_, std::forward<FILTER>(filter));
    }

    /*
     * Performs a rectangular window query. The parameters are the min and max keys which
     * contain the minimum respectively the maximum keys in every dimension.
     * @param query_box The query window.
     * @param filter An optional filter function. The filter function allows filtering entries and
     * sub-nodes before they are returned or traversed. Any filter function must follow the
     * signature of the default 'FilterNoOp`.
     * @return Result iterator.
     */
    template <typename FILTER = FilterNoOp>
    auto begin_query(
        const PhBox<DIM, ScalarInternal>& query_box, FILTER&& filter = FILTER()) const {
        auto pair = find_starting_node(query_box);
        return IteratorHC<T, CONVERT, FILTER>(
            *pair.first,
            query_box.min(),
            query_box.max(),
            converter_,
            std::forward<FILTER>(filter));
    }

    /*
     * Locate nearest neighbors for a given point in space.
     *
     * Example for distance function: auto fn = DistanceEuclidean<DIM>
     * auto iter = tree.begin_knn_query<DistanceEuclidean<DIM>>
     *
     * @param min_results number of entries to be returned. More entries may or may not be returned
     * when several entries have the same distance.
     * @param center center point
     * @param distance_function optional distance function, defaults to euclidean distance
     * @param filter optional filter predicate that excludes nodes/entries before their distance is
     * calculated.
     * @return Result iterator.
     */
    template <typename DISTANCE, typename FILTER = FilterNoOp>
    auto begin_knn_query(
        size_t min_results,
        const KeyT& center,
        DISTANCE&& distance_function = DISTANCE(),
        FILTER&& filter = FILTER()) const {
        return IteratorKnnHS<T, CONVERT, DISTANCE, FILTER>(
            root_,
            min_results,
            center,
            converter_,
            std::forward<DISTANCE>(distance_function),
            std::forward<FILTER>(filter));
    }

    /*
     * @return An iterator representing the tree's 'end'.
     */
    auto end() const {
        return IteratorEnd<EntryT>();
    }

    /*
     * Remove all entries from the tree.
     */
    void clear() {
        num_entries_ = 0;
        root_ = EntryT({}, NodeT{}, detail::MAX_BIT_WIDTH<ScalarInternal> - 1);
    }

    /*
     * @return the number of entries (key/value pairs) in the tree.
     */
    [[nodiscard]] size_t size() const {
        return num_entries_;
    }

    /*
     * @return 'true' if the tree is empty, otherwise 'false'.
     */
    [[nodiscard]] bool empty() const {
        return num_entries_ == 0;
    }

  private:
    /*
     * This function is only for debugging.
     */
    auto GetDebugHelper() const {
        return DebugHelperV20(root_, num_entries_);
    }

    /*
     * Motivation: Point queries a la find() are faster than window queries.
     * Since a window query may have a significant common prefix in their min and max coordinates,
     * the part with the common prefix can be executed as point query.
     *
     * This works if there really is a common prefix, e.g. when querying point data or when
     * querying box data with QueryInclude. Unfortunately, QueryIntersect queries have +/-0 infinity
     * in their coordinates, so their never is an overlap.
     */
    auto find_starting_node(const PhBox<DIM, ScalarInternal>& query_box) const {
        auto& prefix = query_box.min();
        bit_width_t max_conflicting_bits =
            detail::NumberOfDivergingBits(query_box.min(), query_box.max());
        const EntryT* parent = &root_;
        if (max_conflicting_bits > root_.GetNodePostfixLen()) {
            // Abort early if we have no shared prefix in the query
            return std::make_pair(&root_, root_.GetNode().Entries().cend());
        }
        auto entry_iter =
            root_.GetNode().FindPrefix(prefix, max_conflicting_bits, root_.GetNodePostfixLen());
        while (entry_iter != parent->GetNode().Entries().end() && entry_iter->second.IsNode() &&
               entry_iter->second.GetNodePostfixLen() >= max_conflicting_bits) {
            parent = &entry_iter->second;
            entry_iter = parent->GetNode().FindPrefix(
                prefix, max_conflicting_bits, parent->GetNodePostfixLen());
        }
        return std::make_pair(parent, entry_iter);
    }

    size_t num_entries_;
    // Contract: root_ contains a Node with 0 or more entries. The root node is the only Node
    // that is allowed to have less than two entries.
    EntryT root_;
    CONVERT* converter_;
};

}  // namespace improbable::phtree::v20

#endif  // PHTREE_V20_PHTREE_V20_H
