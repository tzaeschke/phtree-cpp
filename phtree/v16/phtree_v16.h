/*
 * Copyright 2020 Improbable Worlds Limited
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

#ifndef PHTREE_V16_PHTREE_V16_H
#define PHTREE_V16_PHTREE_V16_H

#include "debug_helper_v16.h"
#include "for_each.h"
#include "for_each_hc.h"
#include "iterator_full.h"
#include "iterator_hc.h"
#include "iterator_knn_hs.h"
#include "iterator_with_parent.h"
#include "node.h"

namespace improbable::phtree::v16 {

/*
 * The PH-Tree is an ordered index on an n-dimensional space (quad-/oct-/2^n-tree) where each
 * dimension is (by default) indexed by a 64 bit integer. The index ordered follows z-order / Morton
 * order. The index is effectively a 'map', i.e. each key is associated with at most one value.
 *
 * Keys are points in n-dimensional space.
 *
 * This API behaves similar to std::map, see function descriptions for details.
 *
 * Loosely based on PH-Tree Java, V16, see http://www.phtree.org
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
class PhTreeV16 {
    friend PhTreeDebugHelper;
    using ScalarExternal = typename CONVERT::ScalarExternal;
    using ScalarInternal = typename CONVERT::ScalarInternal;
    using KeyT = typename CONVERT::KeyInternal;
    using EntryT = Entry<DIM, T, ScalarInternal>;

  public:
    static_assert(!std::is_reference<T>::value, "Reference type value are not supported.");
    static_assert(std::is_signed<ScalarInternal>::value, "ScalarInternal must be a signed type");
    static_assert(
        std::is_integral<ScalarInternal>::value, "ScalarInternal must be an integral type");
    static_assert(
        std::is_arithmetic<ScalarExternal>::value, "ScalarExternal must be an arithmetic type");
    static_assert(DIM >= 1 && DIM <= 63, "This PH-Tree supports between 1 and 63 dimensions");

    PhTreeV16(CONVERT* converter)
    : num_entries_{0}, root_{MAX_BIT_WIDTH<ScalarInternal> - 1}, converter_{converter} {}

    PhTreeV16(const PhTreeV16& other) = delete;
    PhTreeV16& operator=(const PhTreeV16& other) = delete;
    PhTreeV16(PhTreeV16&& other) noexcept = default;
    PhTreeV16& operator=(PhTreeV16&& other) noexcept = default;
    ~PhTreeV16() noexcept = default;

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
        auto* current_entry = &root_;
        bool is_inserted = false;
        while (current_entry->IsNode()) {
            current_entry = &current_entry->GetNode().Emplace(
                is_inserted, key, current_entry->GetNodePostfixLen(), std::forward<Args>(args)...);
        }
        num_entries_ += is_inserted;
        return {current_entry->GetValue(), is_inserted};
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

            if (!iterator.GetParentNodeEntry()) {
                // No hint available, use standard try_emplace()
                return try_emplace(key, std::forward<Args>(args)...);
            }

            auto* parent_entry = iterator.GetParentNodeEntry();
            if (NumberOfDivergingBits(key, parent_entry->GetKey()) >
                parent_entry->GetNodePostfixLen() + 1) {
                // replace higher up in the tree
                return try_emplace(key, std::forward<Args>(args)...);
            }

            // replace in node
            auto* current_entry = parent_entry;
            bool is_inserted = false;
            while (current_entry->IsNode()) {
                current_entry = &current_entry->GetNode().Emplace(
                    is_inserted,
                    key,
                    current_entry->GetNodePostfixLen(),
                    std::forward<Args>(args)...);
            }
            num_entries_ += is_inserted;
            return {current_entry->GetValue(), is_inserted};
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
        auto* current_entry = &root_;
        while (current_entry && current_entry->IsNode()) {
            current_entry = current_entry->GetNode().Find(key, current_entry->GetNodePostfixLen());
        }
        return current_entry ? 1 : 0;
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
        if (empty()) {
            return IteratorWithParent<T, CONVERT>(nullptr, nullptr, nullptr, converter_);
        }

        const EntryT* current_entry = &root_;
        const EntryT* current_node = nullptr;
        const EntryT* parent_node = nullptr;
        while (current_entry && current_entry->IsNode()) {
            parent_node = current_node;
            current_node = current_entry;
            current_entry = current_entry->GetNode().Find(key, current_entry->GetNodePostfixLen());
        }

        return IteratorWithParent<T, CONVERT>(current_entry, current_node, parent_node, converter_);
    }

    /*
     * See std::map::erase(). Removes any value associated with the provided key.
     *
     * @return '1' if a value was found, otherwise '0'.
     */
    size_t erase(const KeyT& key) {
        auto* current_entry = &root_;
        // We do not pass in the root entry as parent of a node because we do not want the
        // root entry to be modified. The reason is simply that a lot of the code in this class
        // becomes a lot simpler if we can assume the root entry to contain a node.
        EntryT* non_root_current_entry = nullptr;
        bool found = false;
        while (current_entry) {
            auto* child_entry = current_entry->GetNode().Erase(
                key, non_root_current_entry, current_entry->GetNodePostfixLen(), found);
            current_entry = child_entry;
            non_root_current_entry = child_entry;
        }
        num_entries_ -= found;
        return found;
    }

    /*
     * See std::map::erase(). Removes any value at the given iterator location.
     *
     *
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
        if constexpr (std::is_same_v<ITERATOR, IteratorWithParent<T, CONVERT>>) {
            const auto& iter_rich = static_cast<const IteratorWithParent<T, CONVERT>&>(iterator);
            if (!iter_rich.GetCurrentNodeEntry() || iter_rich.GetCurrentNodeEntry() == &root_) {
                // Do _not_ use the root entry, see erase(key). Start searching from the top.
                return erase(iter_rich.GetCurrentResult()->GetKey());
            }
            bool found = false;
            assert(iter_rich.GetCurrentNodeEntry() && iter_rich.GetCurrentNodeEntry()->IsNode());
            iter_rich.GetCurrentNodeEntry()->GetNode().Erase(
                iter_rich.GetCurrentResult()->GetKey(),
                iter_rich.GetCurrentNodeEntry(),
                iter_rich.GetCurrentNodeEntry()->GetNodePostfixLen(),
                found);
            num_entries_ -= found;
            return found;
        }
        // There may be no entry because not every iterator sets it.
        return erase(iterator.GetCurrentResult()->GetKey());
    }

    /*
     * Iterates over all entries in the tree. The optional filter allows filtering entries and nodes
     * (=sub-trees) before returning / traversing them. By default all entries are returned. Filter
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
        ForEachHC<T, CONVERT, CALLBACK, FILTER>(
            query_box.min(),
            query_box.max(),
            converter_,
            std::forward<CALLBACK>(callback),
            std::forward<FILTER>(filter))
            .Traverse(root_);
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
        return IteratorHC<T, CONVERT, FILTER>(
            root_, query_box.min(), query_box.max(), converter_, std::forward<FILTER>(filter));
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
        root_ = EntryT(MAX_BIT_WIDTH<ScalarInternal> - 1);
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
        return DebugHelperV16(root_, num_entries_);
    }

  private:
    size_t num_entries_;
    // Contract: root_ contains a Node with 0 or more entries (the root node is the only Node
    // that is allowed to have less than two entries.
    EntryT root_;
    CONVERT* converter_;
};

}  // namespace improbable::phtree::v16

#endif  // PHTREE_V16_PHTREE_V16_H
