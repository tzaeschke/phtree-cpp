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

#ifndef PHTREE_PHTREE_GRID_INDEX_H
#define PHTREE_PHTREE_GRID_INDEX_H

#include "common/b_plus_tree_hash_map.h"
#include "common/common.h"
#include "phtree_multimap.h"
#include <unordered_set>

namespace improbable::phtree {

/*
 * PH-Tree grid-index main class.
 *
 * The PhTreeGridIndex is a wrapper around a normal PH-tree multi-map.
 * The grid-index has much faster relocate() operations: In case of small movements the
 * cost is O(1) (basically O(0)!). The tree can see whether an entry would stay in the same bin,
 * if it does, the tree is not traversed, the cost of the operation is mainly comparing the old
 * and new key plus some maths.
 *
 * Internally, the grid index just rounds the coordinates to a configurable grid. That's it.
 *
 * The API follows mostly the std::unordered_multimap, exceptions are pointed out.
 * Differences to PhTree
 * - This is a multi-map and hence follows the std::unordered_multimap rather than std::map
 * - erase() returns an iterator instead of a pairs {iterator, bool)
 * - similar to the normal PH-Tree, emplace() returns a reference to the value instead of an
 * iterator
 *
 * For more information please refer to the README of this project.
 */

namespace {

template <dimension_t DIM, typename SCALAR_EXTERNAL, typename SCALAR_INTERNAL>
class ConverterGridIndex : public ConverterPointBase<DIM, SCALAR_EXTERNAL, SCALAR_INTERNAL> {
    using BASE = ConverterPointBase<DIM, SCALAR_EXTERNAL, SCALAR_INTERNAL>;

  public:
    using Point = typename BASE::KeyExternal;
    using PointInternal = typename BASE::KeyInternal;
    using QueryBox = typename BASE::QueryBoxExternal;
    using QueryBoxInternal = typename BASE::QueryBoxInternal;
    using ScalarExternal = typename BASE::ScalarExternal;
    using ScalarInternal = typename BASE::ScalarInternal;

  public:
    explicit ConverterGridIndex(double cell_edge_length)
    : post_{cell_edge_length}, pre_{1. / cell_edge_length} {}

    ConverterGridIndex(const ConverterGridIndex& other) = default;
    ConverterGridIndex& operator=(const ConverterGridIndex& other) = default;
    ConverterGridIndex(ConverterGridIndex&& other) noexcept = default;
    ConverterGridIndex& operator=(ConverterGridIndex&& other) noexcept = default;
    ~ConverterGridIndex() noexcept = default;

    [[nodiscard]] PointInternal pre(const Point& point) const {
        PointInternal p{};
        for (dimension_t d = 0; d < DIM; ++d) {
            p[d] = static_cast<ScalarInternal>(point[d] * pre_);
        }
        return p;
    }

    [[nodiscard]] Point post(const PointInternal& in) const {
        Point p{};
        for (dimension_t d = 0; d < DIM; ++d) {
            p[d] = static_cast<ScalarExternal>(in[d] * post_);
        }
        return p;
    }

    [[nodiscard]] QueryBoxInternal pre_query(const QueryBox& box) const {
        return {pre(box.min()), pre(box.max())};
    }

  private:
    const double post_;
    const double pre_;
};

/*
 * Base class for the internal PH-Tree multi-map iterators.
 *
 * This base class must be distinct from the other Iterator classes because it must be agnostic of
 * the types of the fields that hold iterators. If it knew about these types then we would need
 * to provide them for the ==/!= operators, which would then make it impossible to compare
 * the generic end() iterator with any specialized iterator.
 */
// TODO merge with "Normal" ?!?!?
template <typename PHTREE>
class IteratorBaseGI {
    friend PHTREE;
    using T = typename PHTREE::ValueType;

  public:
    explicit IteratorBaseGI() noexcept : current_value_ptr_{nullptr} {}

    T& operator*() const noexcept {
        assert(current_value_ptr_);
        return const_cast<T&>(*current_value_ptr_);
    }

    T* operator->() const noexcept {
        assert(current_value_ptr_);
        return const_cast<T*>(current_value_ptr_);
    }

    friend bool operator==(
        const IteratorBaseGI<PHTREE>& left, const IteratorBaseGI<PHTREE>& right) noexcept {
        return left.current_value_ptr_ == right.current_value_ptr_;
    }

    friend bool operator!=(
        const IteratorBaseGI<PHTREE>& left, const IteratorBaseGI<PHTREE>& right) noexcept {
        return left.current_value_ptr_ != right.current_value_ptr_;
    }

  protected:
    void SetFinished() noexcept {
        current_value_ptr_ = nullptr;
    }

    void SetCurrentValue(const T* current_value_ptr) noexcept {
        current_value_ptr_ = current_value_ptr;
    }

  private:
    const T* current_value_ptr_;
};

template <typename ITERATOR_PH, typename PHTREE, typename FILTER>
class IteratorNormalGI : public IteratorBaseGI<PHTREE> {
    friend PHTREE;

  public:
    explicit IteratorNormalGI() noexcept : IteratorBaseGI<PHTREE>(), iter_ph_{} {}

    template <typename ITER_PH, typename FILT>
    IteratorNormalGI(ITER_PH&& iter_ph, FILT&& filter) noexcept
    : IteratorBaseGI<PHTREE>()
    , iter_ph_{std::forward<ITER_PH>(iter_ph)}
    , filter_{std::forward<FILT>(filter)} {
        FindNextElement();
    }

    IteratorNormalGI& operator++() noexcept {
        ++iter_ph_;
        FindNextElement();
        return *this;
    }

    IteratorNormalGI operator++(int) noexcept {
        IteratorNormalGI iterator(this->iter_ph_, filter_);  // TODO ... ?
        ++(*this);
        return iterator;
    }

    /*
     * Returns the external key (the 'first' part of the key/value pair).
     */
    auto first() const {
        return iter_ph_.first();
    }

  protected:
    auto& GetIteratorOfPhTree() const noexcept {
        return iter_ph_;
    }

  private:
    void FindNextElement() {
        while (!iter_ph_.__is_end()) {
            // We filter only entries here, nodes are filtered elsewhere
            auto& entry = *iter_ph_;
            // TODO filter
            if (filter_(entry.first)) {
                this->SetCurrentValue(&(entry.second));
                return;
            }
            ++iter_ph_;
        }
        // finished
        this->SetFinished();
    }

    ITERATOR_PH iter_ph_;
    FILTER filter_;
};

template <typename ITERATOR_PH, typename PHTREE, typename FILTER>
class IteratorKnnGI : public IteratorNormalGI<ITERATOR_PH, PHTREE, FILTER> {
  public:
    template <typename ITER_PH, typename FILT>
    IteratorKnnGI(ITER_PH&& iter_ph, FILT&& filter) noexcept
    : IteratorNormalGI<ITER_PH, PHTREE, FILTER>(
          std::forward<ITER_PH>(iter_ph), std::forward<FILT>(filter)) {}

    [[nodiscard]] double distance() const noexcept {
        return this->GetIteratorOfPhTree().distance();
    }
};

}  // namespace

template <typename Key, typename T>
using PhTreeGridIndexEntry = std::pair<Key, T>;
}  // namespace improbable::phtree

namespace std {
// template <>
template <typename Key, typename T>
struct hash<typename improbable::phtree::PhTreeGridIndexEntry<Key, T>> {
    size_t operator()(const typename improbable::phtree::PhTreeGridIndexEntry<Key, T>& x) const {
        return std::hash<T>{}(x.second);
    }
};
};  // namespace std

namespace improbable::phtree {
/*
 * The PhTreeMultiMap class.
 */
template <
    dimension_t DIM,
    typename T,
    // typename CONVERTER = ConverterNoOp<DIM, scalar_64_t>,
    typename CONVERTER = ConverterGridIndex<DIM, double, scalar_64_t>,
    typename BUCKET = b_plus_tree_hash_set<T>,
    bool POINT_KEYS = true,
    typename DEFAULT_QUERY_TYPE = QueryPoint>
class PhTreeGridIndex {
    using KeyInternal = typename CONVERTER::KeyInternal;
    using Key = typename CONVERTER::KeyExternal;
    static constexpr dimension_t DimInternal = CONVERTER::DimInternal;
    using PHTREE = PhTreeGridIndex<DIM, T, CONVERTER, BUCKET, POINT_KEYS, DEFAULT_QUERY_TYPE>;
    using ValueType = T;
    using EndType = decltype(std::declval<PhTreeMultiMap<
                                 DIM,
                                 PhTreeGridIndexEntry<Key, T>,
                                 CONVERTER,
                                 BUCKET,
                                 POINT_KEYS,
                                 DEFAULT_QUERY_TYPE>>()
                                 .end());

    friend PhTreeDebugHelper;
    friend IteratorBaseGI<PHTREE>;

  public:
    using QueryBox = typename CONVERTER::QueryBoxExternal;
    using EntryT = PhTreeGridIndexEntry<Key, T>;

  private:
    using BUCKET_Internal = b_plus_tree_hash_set<EntryT>;

  public:
    explicit PhTreeGridIndex(double cell_edge_length = 100) : tree_{CONVERTER{cell_edge_length}} {}

    explicit PhTreeGridIndex(CONVERTER converter) : tree_{converter} {}

    PhTreeGridIndex(const PhTreeGridIndex& other) = delete;
    PhTreeGridIndex& operator=(const PhTreeGridIndex& other) = delete;
    PhTreeGridIndex(PhTreeGridIndex&& other) noexcept = default;
    PhTreeGridIndex& operator=(PhTreeGridIndex&& other) noexcept = default;
    ~PhTreeGridIndex() noexcept = default;

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
        auto result = tree_.try_emplace(key, EntryT{key, std::forward<Args>(args)...});
        return {const_cast<T&>(result.first.second), result.second};
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
        auto result = tree_.emplace_hint(
            iterator.GetIteratorOfPhTree(), key, EntryT{key, std::forward<Args>(args)...});
        return {const_cast<T&>(result.first.second), result.second};
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
        auto iter = tree_.find(key);
        size_t n = 0;
        while (iter != tree_.end()) {
            n += (key == iter->first);
            ++iter;
        }
        return n;
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
        return tree_.estimate_count(query_box, query_type);
        //        // TODO... use box filter
        //        size_t n = 0;
        //        auto counter_lambda = [&](const Key&, const EntryT& bucket) { ++n; };
        //        // auto filter = [&](const Key&, const BUCKET& bucket) { n += bucket.size(); };
        //        tree_.for_each(query_box, counter_lambda, FilterNoOp{}, query_type);
        //        return n;
    }

    /*
     * See std::unordered_multimap::find().
     *
     * @param key the key to look up
     * @return an iterator that points either to the first value associated with the key or
     * to {@code end()} if no value was found
     */
    auto find(const Key& key) const {
        auto filter = [key = key](const Key& key2) noexcept { return key == key2; };
        return CreateIterator(tree_.find(key), std::move(filter));
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
        auto filter = [key = key](const Key& key2) noexcept { return key == key2; };
        return CreateIterator(tree_.find(key, create(key, value)), std::move(filter));
    }

    /*
     * See std::unordered_multimap::erase(). Removes the provided key/value pair if it exists.
     *
     * @return '1' if the key/value pair was found, otherwise '0'.
     */
    size_t erase(const Key& key, const T& value) {
        return tree_.erase(key, create(key, value));
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
        return tree_.erase(iterator.GetIteratorOfPhTree());
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
    size_t relocate(
        const Key& old_key, const Key& new_key, T2&& value, bool verify_exists = false) {
        // TODO document verify_exists,
        // TODO do we need to check coordinates? Document this!!
        // TODO update old/new key? With verify=false we can ignore updating the key!!
        // return tree_.relocate(old_key, new_key, create(old_key, value), verify_exists);

        auto update_fn = [&value, &old_key, &new_key](const EntryT& e) -> size_t {
            if (e.second == value && e.first == old_key) {
                const_cast<Key&>(e.first) = new_key;
                return true;
            }
            return false;
        };
        return tree_.relocate_if(old_key, new_key, std::move(update_fn), true);

        //        auto fn = [&value](BUCKET& src, BUCKET& dst) -> size_t {
        //            auto it = src.find(value);
        //            if (it != src.end() && dst.emplace(std::move(*it)).second) {
        //                src.erase(it);
        //                return 1;
        //            }
        //            return 0;
        //        };
        //        auto count_fn = [&value](BUCKET& src) -> size_t { return src.find(value) !=
        //        src.end(); }; return tree_._relocate_mm(
        //            converter_.pre(old_key), converter_.pre(new_key), verify_exists, fn,
        //            count_fn);
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
     *              ConverterMultiply.
     *
     * @return the number of values that were relocated.
     */
    template <typename PREDICATE>
    size_t relocate_if(
        const Key& old_key, const Key& new_key, PREDICATE&& pred_fn, bool verify_exists = false) {
        // TODO document verify_exists,
        // TODO do we need to check coordinates? Document this!!
        return tree_.relocate_if(old_key, new_key, std::forward<PREDICATE>(pred_fn), verify_exists);
        //        auto fn = [&pred_fn](BUCKET& src, BUCKET& dst) -> size_t {
        //            size_t result = 0;
        //            auto iter_src = src.begin();
        //            while (iter_src != src.end()) {
        //                if (pred_fn(*iter_src) && dst.emplace(std::move(*iter_src)).second) {
        //                    iter_src = src.erase(iter_src);
        //                    ++result;
        //                } else {
        //                    ++iter_src;
        //                }
        //            }
        //            return result;
        //        };
        //        auto count_fn = [&pred_fn](BUCKET& src) -> size_t {
        //            size_t result = 0;
        //            auto iter_src = src.begin();
        //            while (iter_src != src.end()) {
        //                if (pred_fn(*iter_src)) {
        //                    ++result;
        //                }
        //                ++iter_src;
        //            }
        //            return result;
        //        };
        //        return tree_._relocate_mm(
        //            converter_.pre(old_key), converter_.pre(new_key), verify_exists, fn,
        //            count_fn);
    }

    /*
     * Relocates all values from one coordinate to another.
     * Returns an iterator pointing to the relocated data (or end(), if the relocation failed).
     */
    auto relocate_all(const Key& old_key, const Key& new_key) {
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
        tree_.for_each(
            NoOpCallback{},
            WrapCallbackFilter<CALLBACK, FILTER>{
                std::forward<CALLBACK>(callback), std::forward<FILTER>(filter), converter()});
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
        tree_.template for_each<NoOpCallback, WrapCallbackFilterQuery<CALLBACK, FILTER>>(
            query_box,
            {},
            {std::forward<CALLBACK>(callback),
             std::forward<FILTER>(filter),
             converter(),
             query_box},
            query_type);
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
        return CreateIterator(tree_.begin(WrapFilter<FILTER>(std::forward<FILTER>(filter))));
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
        auto key_filter = [query_box = query_box](const Key& key) noexcept {
            auto& min = query_box.min();
            auto& max = query_box.max();
            for (dimension_t d = 0; d < DIM; ++d) {
                if (key[d] < min[d] || key[d] > max[d]) {
                    return false;
                }
            }
            return true;
        };
        return CreateIterator(
            tree_.begin_query(
                query_box, WrapFilter<FILTER>(std::forward<FILTER>(filter)), query_type),
            std::move(key_filter));
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
        // TODO filter
        return CreateIteratorKnn(tree_.begin_knn_query(
            min_results,
            center,
            std::forward<DISTANCE>(distance_function),
            WrapFilter<FILTER>(std::forward<FILTER>(filter))));
    }

    /*
     * @return An iterator representing the tree's 'end'.
     */
    auto end() const {
        return IteratorNormalGI<EndType, PHTREE, NoOpFilterGI>{};
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
        return tree_.converter();
    }

  private:
    // This is used by PhTreeDebugHelper
    const auto& GetInternalTree() const {
        return tree_.GetInternalTree();
    }

    void CheckConsistencyExternal() const {
        tree_.CheckConsistencyExternal();
    }

    template <typename... Args>
    EntryT create(const Key& key, Args&&... args) const {
        return std::make_pair(key, std::forward<Args>(args)...);
    }

    struct NoOpCallback {
        constexpr void operator()(const Key&, const EntryT&) const noexcept {}
    };

    struct NoOpFilterGI {
        constexpr bool operator()(const Key&) const noexcept {
            return true;
        }
    };

    template <typename OUTER_ITER, typename KEY_FILTER = NoOpFilterGI>
    auto CreateIterator(OUTER_ITER&& outer_iter, KEY_FILTER&& filter = KEY_FILTER()) const {
        return IteratorNormalGI<OUTER_ITER, PHTREE, KEY_FILTER>(
            std::forward<OUTER_ITER>(outer_iter), std::forward<KEY_FILTER>(filter));
    }

    template <typename OUTER_ITER, typename KEY_FILTER = NoOpFilterGI>
    auto CreateIteratorKnn(OUTER_ITER&& outer_iter, KEY_FILTER&& filter = KEY_FILTER()) const {
        return IteratorKnnGI<OUTER_ITER, PHTREE, KEY_FILTER>(
            std::forward<OUTER_ITER>(outer_iter), std::forward<KEY_FILTER>(filter));
    }

    template <typename FILTER>
    class WrapFilter {
      public:
        template <typename F>
        WrapFilter(F&& filter) : filter_{std::forward<F>(filter)} {}

        template <typename BucketT>
        [[nodiscard]] constexpr bool IsEntryValid(const KeyInternal&, const BucketT& e) const {
            return true;  // filter_.IsEntryValid(e.first, e.second);
        }
        [[nodiscard]] constexpr bool IsNodeValid(const KeyInternal&, int) const {
            // TODO? Remove filter methods for grid ?!?!
            return true;
        }
        [[nodiscard]] constexpr bool IsBucketEntryValid(
            const KeyInternal& k, const EntryT& e) const {
            // TODO avoid using key-internal
            return filter_.IsBucketEntryValid(k, e.second);
        }

      private:
        FILTER filter_;
    };

    /*
     * This wrapper wraps the Filter and Callback such that the callback is called for every
     * entry in any bucket that matches the user defined IsEntryValid().
     */
    template <typename CALLBACK, typename FILTER>
    class WrapCallbackFilter {
      public:
        /*
         * We always have two iterators, one that traverses the PH-Tree and returns 'buckets', the
         * other iterator traverses the returned buckets.
         * The wrapper ensures that the callback is called for every entry in a bucket..
         */
        template <typename CB, typename F>
        WrapCallbackFilter(CB&& callback, F&& filter, const CONVERTER& converter)
        : callback_{std::forward<CB>(callback)}
        , filter_{std::forward<F>(filter)}
        , converter_{converter} {}

        [[nodiscard]] inline bool IsEntryValid(
            const KeyInternal& internal_key, const BUCKET& bucket) {
            // TODO???
            //   We can roughly filter the bucket by key, but we need to traverse all
            //   entries anyway to get the correct key.
            //   Problem: we cannot easily map the type of the internal bucket to the external
            //   bucket because of the different Entry type.
            //   However we can simply forward the modified bucket type, it is easy to use,
            //   even if it does not comply with normal signature....

            //            if (filter_.IsEntryValid(internal_key, bucket)) {
            //                auto key = converter_.post(internal_key);
            //                for (auto& entry : bucket) {
            //                    if (filter_.IsBucketEntryValid(internal_key, entry)) {
            //                        callback_(key, entry);
            //                    }
            //                }
            //            }
            //            // Return false. We already called the callback.
            //            return false;
            return true;
        }

        template <typename ValueT>
        [[nodiscard]] inline bool IsBucketEntryValid(
            const KeyInternal& internal_key, const ValueT& entry) const noexcept {
            if (filter_.IsBucketEntryValid(internal_key, entry.second)) {
                callback_(entry.first, entry.second);
                return true;
            }
            // Return false. We already called the callback.
            return false;
        }

        [[nodiscard]] inline bool IsNodeValid(const KeyInternal& prefix, int bits_to_ignore) {
            // TODO document this?!? We cannot check the nodes.....
            // TODO disable filters all together?
            return true;
            // return filter_.IsNodeValid(prefix, bits_to_ignore);
        }

      private:
        CALLBACK callback_;
        FILTER filter_;
        const CONVERTER& converter_;
    };

    template <typename CALLBACK, typename FILTER>
    class WrapCallbackFilterQuery {
      public:
        /*
         * We always have two iterators, one that traverses the PH-Tree and returns 'buckets', the
         * other iterator traverses the returned buckets.
         * The wrapper ensures that the callback is called for every entry in a bucket..
         */
        template <typename CB, typename F>
        WrapCallbackFilterQuery(
            CB&& callback, F&& filter, const CONVERTER& converter, const QueryBox& query)
        : callback_{std::forward<CB>(callback)}
        , filter_{std::forward<F>(filter)}
        , converter_{converter}
        , query_{query} {}

        [[nodiscard]] inline bool IsEntryValid(const KeyInternal&, const BUCKET_Internal&) {
            return true;
        }

        template <typename ValueT>
        [[nodiscard]] inline bool IsBucketEntryValid(
            const KeyInternal& internal_key, const ValueT& entry) const noexcept {
            auto& min = query_.min();
            auto& max = query_.max();
            auto& key = entry.first;
            for (dimension_t d = 0; d < DIM; ++d) {
                if (key[d] < min[d] || key[d] > max[d]) {
                    return false;
                }
            }

            if (filter_.IsBucketEntryValid(internal_key, entry.second)) {
                callback_(entry.first, entry.second);
                return true;
            }
            // Return false. We already called the callback.
            return false;
        }

        [[nodiscard]] inline bool IsNodeValid(const KeyInternal& prefix, int bits_to_ignore) {
            // TODO document this?!? We cannot check the nodes.....
            // TODO disable filters all together?
            return true;
            // return filter_.IsNodeValid(prefix, bits_to_ignore);
        }

      private:
        CALLBACK callback_;
        FILTER filter_;
        const CONVERTER& converter_;
        QueryBox query_;
    };

    PhTreeMultiMap<DIM, EntryT, CONVERTER, BUCKET_Internal, POINT_KEYS, DEFAULT_QUERY_TYPE> tree_;
};

/**
 * A PH-Tree multi-map that uses (axis aligned) points as keys.
 * The points are defined with 64bit 'double' floating point coordinates.
 *
 * See 'PhTreeD' for details.
 */
template <
    dimension_t DIM,
    typename T,
    typename CONVERTER = ConverterGridIndex<DIM, double, scalar_64_t>,
    // TODO !!!!!!!!!!!!!!!!11
    typename BUCKET = b_plus_tree_hash_set<PhTreeGridIndexEntry<PhPointD<DIM>, T>>>
using PhTreeGridIndexD = PhTreeGridIndex<DIM, T, CONVERTER, BUCKET>;

template <
    dimension_t DIM,
    typename T,
    typename CONVERTER_BOX,
    typename BUCKET = b_plus_tree_hash_set<PhTreeGridIndexEntry<PhPointD<DIM>, T>>>
using PhTreeGridIndexBox = PhTreeGridIndex<DIM, T, CONVERTER_BOX, BUCKET, false, QueryIntersect>;

/**
 * A PH-Tree multi-map that uses (axis aligned) boxes as keys.
 * The boxes are defined with 64bit 'double' floating point coordinates.
 *
 * See 'PhTreeD' for details.
 */
template <
    dimension_t DIM,
    typename T,
    typename CONVERTER_BOX = ConverterBoxIEEE<DIM>,
    typename BUCKET = b_plus_tree_hash_set<PhTreeGridIndexEntry<PhPointD<DIM>, T>>>
using PhTreeGridIndexBoxD = PhTreeGridIndexBox<DIM, T, CONVERTER_BOX, BUCKET>;

}  // namespace improbable::phtree

#endif  // PHTREE_PHTREE_GRID_INDEX_H
