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

#ifndef PHTREE_PHTREE_MULTIMAP2_H
#define PHTREE_PHTREE_MULTIMAP2_H

#include "common/b_plus_tree_hash_map.h"
#include "common/common.h"
#include "v20/phtree_v20.h"
#include <unordered_set>

namespace improbable::phtree {

/*
 * PH-Tree multi-map main class.
 *
 * The PhTreeMultiMap is a wrapper around a normal PH-Tree (single value per key). The wrapper uses
 * collections to store more than one value per key.
 * By default, this multi-map is backed by std::unordered_set<T>.
 *
 * The API follows mostly the std::unordered_multimap, exceptions are pointed out.
 *
 * Differences to PhTree
 * - This is a multi-map and hence follows the std::unordered_multimap rather than std::map
 * - erase() returns an iterator instead of a pairs {iterator, bool)
 * - similar to the normal PH-Tree, emplace() returns a reference to the value instead of an
 * iterator
 *
 * For more information please refer to the README of this project.
 */

/*
 * The PhTreeMultiMap class.
 */
template <
    dimension_t DIM,
    typename T,
    typename CONVERTER = ConverterNoOp<DIM, scalar_64_t>,
    bool POINT_KEYS = true,
    typename DEFAULT_QUERY_TYPE = QueryPoint>
class PhTreeMultiMap2 {
    using ScalarInternal = typename CONVERTER::ScalarInternal;
    using KeyInternal = typename CONVERTER::KeyInternal;
    using Key = typename CONVERTER::KeyExternal;
    static constexpr dimension_t DimInternal = CONVERTER::DimInternal;
    using PHTREE = PhTreeMultiMap2<DIM, T, CONVERTER, POINT_KEYS, DEFAULT_QUERY_TYPE>;
    using ValueType = T;
    // using IterType = decltype(std::declval<T>().begin());
    //  TODO use auto instead...?
    using IterType = decltype(std::declval<v20::PhTreeV20<DimInternal, T, CONVERTER>>().begin());
    using EndType = decltype(std::declval<v20::PhTreeV20<DimInternal, T, CONVERTER>>().end());

    friend PhTreeDebugHelper;

  public:
    using QueryBox = typename CONVERTER::QueryBoxExternal;

    explicit PhTreeMultiMap2(CONVERTER converter = CONVERTER())
    : tree_{&converter_}, converter_{converter} {}

    PhTreeMultiMap2(const PhTreeMultiMap2& other) = delete;
    PhTreeMultiMap2& operator=(const PhTreeMultiMap2& other) = delete;
    PhTreeMultiMap2(PhTreeMultiMap2&& other) noexcept = default;
    PhTreeMultiMap2& operator=(PhTreeMultiMap2&& other) noexcept = default;
    ~PhTreeMultiMap2() noexcept = default;

    /*
     *  Attempts to build and insert a key and a value into the tree.
     *
     *  @param key The key for the new entry.
     *
     *  @param args  Arguments used to generate a new value.
     *
     *  @return  A pair, whose first element points to the possibly inserted pair,
     *           and whose second element is a bool that is true if the pair was actually inserted.
     *
     * This function attempts to build and insert a (key, value) pair into the tree. The PH-Tree is
     * effectively a multi-set, so if an entry with the same key/value was already in the tree, it
     * returns that entry instead of inserting a new one.
     */
    template <typename... Args>
    std::pair<T&, bool> emplace(const Key& key, Args&&... args) {
        auto iter = tree_.try_emplace(converter_.pre(key), std::forward<Args>(args)...);
        //      size_ += iter.second ? 1 : 0;
        return std::move(iter);  // TODO remove...? {const_cast<T&>(*iter.first), iter.second};
    }

    /*
     * The emplace_hint() method uses an iterator as hint for insertion.
     * The hint is ignored if it is not useful or is equal to end().
     *
     * Iterators should normally not be used after the tree has been modified. As an exception to
     * this rule, an iterator can be used as hint if it was previously used with at most one call
     * to erase() and if no other modifications occurred.
     * The following is valid:
     *
     * // Move value from key1 to key2 (if you don't want to use relocate() ).
     * auto iter = tree.find(key1);
     * auto value = iter.second(); // The value may become invalid in erase()
     * erase(iter);
     * emplace_hint(iter, key2, value);  // the iterator can still be used as hint here
     */
    template <typename ITERATOR, typename... Args>
    std::pair<T&, bool> emplace_hint(const ITERATOR& iterator, const Key& key, Args&&... args) {
        auto result = tree_.try_emplace(iterator, converter_.pre(key), std::forward<Args>(args)...);
        // TODO return iterator i.o. pair! See:
        //    https://en.cppreference.com/w/cpp/container/unordered_multimap/emplace
        return result;
    }

    /*
     * See std::unordered_multimap::insert().
     *
     * @return a pair consisting of the inserted value (or to the value that prevented the
     * insertion if the key/value already existed) and a bool denoting whether the insertion
     * took place.
     */
    std::pair<T&, bool> insert(const Key& key, const T& value) {
        return emplace(key, value);
    }

    /*
     * See emplace().
     */
    template <typename... Args>
    std::pair<T&, bool> try_emplace(const Key& key, Args&&... args) {
        return emplace(key, std::forward<Args>(args)...);
    }

    /*
     * See emplace_hint().
     */
    template <typename ITERATOR, typename... Args>
    std::pair<T&, bool> try_emplace(const ITERATOR& iterator, const Key& key, Args&&... args) {
        return emplace_hint(iterator, key, std::forward<Args>(args)...);
    }

    /*
     * @return '1', if a value is associated with the provided key, otherwise '0'.
     */
    size_t count(const Key& key) const {
        return tree_.count(converter_.pre(key));
    }

    /*
     * Estimates the result count of a rectangular window query by counting the sizes of all buckets
     * that overlap with the query box. This estimate function should be much faster than a normal
     * query, especially in trees with many entries per bucket.
     *
     * @param query_box The query window.
     * @param query_type The type of query, such as QueryIntersect or QueryInclude
     */
    template <typename QUERY_TYPE = DEFAULT_QUERY_TYPE>
    size_t estimate_count(QueryBox query_box, QUERY_TYPE query_type = QUERY_TYPE()) const {
        size_t n = 0;
        auto counter_lambda = [&](const Key&, const T& bucket) { ++n; };
        tree_.for_each(query_type(converter_.pre_query(query_box)), counter_lambda);
        return n;
    }

    /*
     * See std::unordered_multimap::find().
     *
     * @param key the key to look up
     * @return an iterator that points either to the first value associated with the key or
     * to {@code end()} if no value was found
     */
    auto find(const Key& key) const {
        return tree_.find(converter_.pre(key));
    }

    /*
     * See std::unordered_multimap::find().
     *
     * @param key the key to look up
     * @param value the value to look up
     * @return an iterator that points either to the associated value of the key/value pair
     * or to {@code end()} if the key/value pair was found
     */
    auto find(const Key& key, const T& value) const {
        return tree_.find(converter_.pre(key), value);
    }

    /*
     * See std::unordered_multimap::erase(). Removes the provided key/value pair if it exists.
     * // TODO document: removes exactly one occurrence!
     *
     * @return '1' if the key/value pair was found, otherwise '0'.
     */
    // TODO do we need this method? It just complicates the API...  Provide erase_if() instead?
    //    -> Not present in std::(unordered_)multimap
    size_t erase(const Key& key, const T& value) {
        return tree_.erase(converter_.pre(key), value);
        //        auto key2 = converter_.pre(key);
        //        auto iter = tree_.find(key2);
        //        size_t n = 0;
        //        // TODO skip first comparison?
        //        while (iter != tree_.end() && iter.GetEntry()->GetKey() == key2) {
        //            if (iter.GetEntry()->GetValue() == value) {
        //                tree_.erase(iter);
        //                --size_;
        //                ++n;
        //            }
        //            ++iter;
        //        }
        //        return n;
    }

    /*
     * See std::map::erase(). Removes any entry located at the provided iterator.
     *
     * This function uses the iterator to directly erase the entry, so it is usually faster than
     * erase(key, value).
     *
     * @return '1' if a value was found, otherwise '0'.
     */
    template <typename ITERATOR>
    size_t erase(const ITERATOR& iterator) {
        // TODO backport to v20
        // std::is_convertible_v<ITERATOR*, IterType*>,
        using IterT1 = decltype(std::declval<v20::PhTreeV20<DimInternal, T, CONVERTER>>().find({}));
        using IterT2 = decltype(std::declval<v20::PhTreeV20<DimInternal, T, CONVERTER>>().end());
        static_assert(
            // std::is_convertible_v<ITERATOR*, v20::IteratorBase<v20::Entry<DIM, T,
            // ScalarInternal>*>>,
            std::is_convertible_v<ITERATOR*, IterT1*> || std::is_convertible_v<ITERATOR*, IterT2*>,
            "erase(iterator) requires an iterator argument. For erasing by key please use "
            "erase(key, value).");
        if (iterator != end()) {
            //            auto& bucket = const_cast<T&>(*iterator.GetIteratorOfPhTree());
            //            size_t old_size = bucket.size();
            //            bucket.erase(iterator.GetIteratorOfBucket());
            //            bool success = bucket.size() < old_size;
            //            if (bucket.empty()) {
            //                success &= tree_.erase(iterator.GetIteratorOfPhTree()) > 0;
            //            }
            //            size_ -= success;
            //            return success;
            return tree_.erase(iterator);
        }
        return 0;
    }

    /*
     * This function attempts to remove the 'value' from 'old_key' and reinsert it for 'new_key'.
     *
     * The relocate function will report _success_ in the following cases:
     * - the value was removed from the old position and reinserted at the new position
     * - the old position and new position are identical.
     *
     * The relocate function will report _failure_ in the following cases:
     * - The value was already present in the new position
     * - The value was not present in the old position
     *
     * In case of _failure_, this function guarantees that the tree remains unchanged
     * or is returned to its original state (i.e. before the function was called).
     *
     * @param old_key The old position
     * @param new_key The new position
     * @param value The value that needs to be relocated. The relocate() method used the value's
     *              '==' operator to identify the entry that should be moved.
     * @param verify_exists This setting toggles whether a relocate() between two identical keys
     *              should verify whether the key actually exist before return '1'.
     *              If set to 'false', this function will return '1' if the keys are identical,
     *              without checking whether the keys actually exist. Avoiding this check can
     *              considerably speed up relocate() calls, especially when using a
     *              ConverterMultiply.
     *
     * @return '1' if a value was found and reinserted, otherwise '0'.
     */
    template <typename T2>
    size_t relocate(const Key& old_key, const Key& new_key, T2&& value, bool verify_exists = true) {
        auto pred_fn = [&value](T& src) -> size_t { return src == value; };
        return tree_.relocate_if(converter_.pre(old_key), converter_.pre(new_key), pred_fn);
    }

    /*
     * This function attempts to remove the 'value' from 'old_key' and reinsert it for 'new_key'.
     *
     * The relocate function will report _success_ in the following cases:
     * - the value was removed from the old position and reinserted at the new position
     * - the old position and new position are identical.
     *
     * The relocate function will report _failure_ in the following cases:
     * - The value was already present in the new position
     * - The value was not present in the old position
     *
     * In case of _failure_, this function guarantees that the tree remains unchanged
     * or is returned to its original state (i.e. before the function was called).
     *
     * @param old_key The old position
     * @param new_key The new position
     * @param predicate The predicate that is used for every value at position old_key to evaluate
     *             whether it should be relocated to new_key.
     * @param verify_exists This setting toggles whether a relocate() between two identical keys
     *              should verify whether the key actually exist before return '1'.
     *              If set to 'false', this function will return '1' if the keys are identical,
     *              without checking whether the keys actually exist. Avoiding this check can
     *              considerably speed up relocate() calls, especially when using a
     *              ConverterMultiply. // TODO remove?
     *
     * @return the number of values that were relocated.
     */
    template <typename PREDICATE>
    size_t relocate_if(
        const Key& old_key, const Key& new_key, PREDICATE&& pred_fn, bool verify_exists = true) {
        return tree_.relocate_if(converter_.pre(old_key), converter_.pre(new_key), pred_fn);
    }

    /*
     * Relocates all values from one coordinate to another.
     * Returns an iterator pointing to the relocated data (or end(), if the relocation failed).
     */
    auto relocate_all(const Key& old_key, const Key& new_key) {
        // TODO remove?
        return tree_.relocate(old_key, new_key);
    }

    /*
     * Iterates over all entries in the tree. The optional filter allows filtering entries and nodes
     * (=sub-trees) before returning / traversing them. By default, all entries are returned. Filter
     * functions must implement the same signature as the default 'FilterNoOp'.
     *
     * @param callback The callback function to be called for every entry that matches the filter.
     * The callback requires the following signature: callback(const PhPointD<DIM> &, const T &)
     * @param filter An optional filter function. The filter function allows filtering entries and
     * sub-nodes before they are passed to the callback or traversed. Any filter function must
     * follow the signature of the default 'FilterNoOp`.
     * The default 'FilterNoOp` filter matches all entries.
     */
    template <typename CALLBACK, typename FILTER = FilterNoOp>
    void for_each(CALLBACK&& callback, FILTER&& filter = FILTER()) const {
        tree_.for_each(std::forward<CALLBACK>(callback), std::forward<FILTER>(filter));
    }

    /*
     * Performs a rectangular window query. The parameters are the min and max keys which
     * contain the minimum respectively the maximum keys in every dimension.
     * @param query_box The query window.
     * @param callback The callback function to be called for every entry that matches the query
     * and filter.
     * The callback requires the following signature: callback(const PhPointD<DIM> &, const T &)
     * @param query_type The type of query, such as QueryIntersect or QueryInclude
     * @param filter An optional filter function. The filter function allows filtering entries and
     * sub-nodes before they are returned or traversed. Any filter function must follow the
     * signature of the default 'FilterNoOp`.
     * The default 'FilterNoOp` filter matches all entries.
     */
    template <
        typename CALLBACK,
        typename FILTER = FilterNoOp,
        typename QUERY_TYPE = DEFAULT_QUERY_TYPE>
    void for_each(
        QueryBox query_box,
        CALLBACK&& callback,
        FILTER&& filter = FILTER(),
        QUERY_TYPE query_type = QUERY_TYPE()) const {
        tree_.for_each(
            query_type(converter_.pre_query(query_box)),
            std::forward<CALLBACK>(callback),
            std::forward<FILTER>(filter));
    }

    /*
     * Iterates over all entries in the tree. The optional filter allows filtering entries and nodes
     * (=sub-trees) before returning / traversing them. By default, all entries are returned. Filter
     * functions must implement the same signature as the default 'FilterNoOp'.
     *
     * @return an iterator over all (filtered) entries in the tree,
     */
    template <typename FILTER = FilterNoOp>
    auto begin(FILTER&& filter = FILTER()) const {
        return tree_.begin(std::forward<FILTER>(filter));
    }

    /*
     * Performs a rectangular window query. The parameters are the min and max keys which
     * contain the minimum respectively the maximum keys in every dimension.
     * @param query_box The query window.
     * @param query_type The type of query, such as QueryIntersect or QueryInclude
     * @param filter An optional filter function. The filter function allows filtering entries and
     * sub-nodes before they are returned or traversed. Any filter function must follow the
     * signature of the default 'FilterNoOp`.
     * @return Result iterator.
     */
    template <typename FILTER = FilterNoOp, typename QUERY_TYPE = DEFAULT_QUERY_TYPE>
    auto begin_query(
        const QueryBox& query_box,
        FILTER&& filter = FILTER(),
        QUERY_TYPE&& query_type = QUERY_TYPE()) const {
        return tree_.begin_query(
            query_type(converter_.pre_query(query_box)), std::forward<FILTER>(filter));
    }

    /*
     * Locate nearest neighbors for a given point in space.
     *
     * NOTE: This method is not (currently) available for box keys.
     *
     * @param min_results number of entries to be returned. More entries may or may not be returned
     * when several entries have the same distance.
     * @param center center point
     * @param distance_function optional distance function, defaults to euclidean distance
     * @param filter optional filter predicate that excludes nodes/entries before their distance is
     * calculated.
     * @return Result iterator.
     */
    template <
        typename DISTANCE,
        typename FILTER = FilterNoOp,
        // Some magic to disable this in case of box keys
        bool DUMMY = POINT_KEYS,
        typename std::enable_if<DUMMY, int>::type = 0>
    auto begin_knn_query(
        size_t min_results,
        const Key& center,
        DISTANCE&& distance_function = DISTANCE(),
        FILTER&& filter = FILTER()) const {
        // We use pre() instead of pre_query() here because, strictly speaking, we want to
        // find the nearest neighbors of a (fictional) key, which may as well be a box.
        return tree_.begin_knn_query(
            min_results,
            converter_.pre(center),
            std::forward<DISTANCE>(distance_function),
            std::forward<FILTER>(filter));
    }

    /*
     * @return An iterator representing the tree's 'end'.
     */
    auto end() const {
        return tree_.end();
    }

    /*
     * Remove all entries from the tree.
     */
    void clear() {
        tree_.clear();
    }

    /*
     * @return the number of entries (key/value pairs) in the tree.
     */
    [[nodiscard]] size_t size() const {
        return tree_.size();
    }

    /*
     * @return 'true' if the tree is empty, otherwise 'false'.
     */
    [[nodiscard]] bool empty() const {
        return tree_.empty();
    }

    /*
     * @return the converter associated with this tree.
     */
    [[nodiscard]] const CONVERTER& converter() const {
        return converter_;
    }

  private:
    // This is used by PhTreeDebugHelper
    const auto& GetInternalTree() const {
        return tree_;
    }

    void CheckConsistencyExternal() const {
        [[maybe_unused]] size_t n = 0;
        for ([[maybe_unused]] const auto& entry : tree_) {
            ++n;
        }
        assert(n == size());
    }

    v20::PhTreeV20<DimInternal, T, CONVERTER> tree_;
    CONVERTER converter_;
};

template <dimension_t DIM, typename T, typename CONVERTER_BOX>
using PhTreeMultiMap2Box = PhTreeMultiMap2<DIM, T, CONVERTER_BOX, false, QueryIntersect>;

/**
 * A PH-Tree multi-map that uses (axis aligned) points as keys.
 * The points are defined with 64bit 'double' floating point coordinates.
 *
 * See 'PhTreeD' for details.
 */
template <dimension_t DIM, typename T, typename CONVERTER = ConverterIEEE<DIM>>
using PhTreeMultiMap2D = PhTreeMultiMap2<DIM, T, CONVERTER>;

/**
 * A PH-Tree multi-map that uses (axis aligned) boxes as keys.
 * The boxes are defined with 64bit 'double' floating point coordinates.
 *
 * See 'PhTreeD' for details.
 */
template <dimension_t DIM, typename T, typename CONVERTER_BOX = ConverterBoxIEEE<DIM>>
using PhTreeMultiMap2BoxD = PhTreeMultiMap2Box<DIM, T, CONVERTER_BOX>;

/**
 * A PH-Tree multi-map that uses (axis aligned) points as keys.
 * The points are defined with 32bit 'float' floating point coordinates.
 *
 * See 'PhTreeF' for details.
 */
template <dimension_t DIM, typename T, typename CONVERTER = ConverterFloatIEEE<DIM>>
using PhTreeMultiMap2F = PhTreeMultiMap2<DIM, T, CONVERTER>;

/**
 * A PH-Tree multi-map that uses (axis aligned) boxes as keys.
 * The boxes are defined with 32bit 'float' floating point coordinates.
 *
 * See 'PhTreeF' for details.
 */
template <dimension_t DIM, typename T, typename CONVERTER_BOX = ConverterBoxFloatIEEE<DIM>>
using PhTreeMultiMapBox2F = PhTreeMultiMap2Box<DIM, T, CONVERTER_BOX>;

}  // namespace improbable::phtree

#endif  // PHTREE_PHTREE_MULTIMAP_H
