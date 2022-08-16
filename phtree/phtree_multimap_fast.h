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

#ifndef PHTREE_PHTREE_MULTIMAP_FAST_H
#define PHTREE_PHTREE_MULTIMAP_FAST_H

#include "common/b_plus_tree_hash_map.h"
#include "common/common.h"
#include "phtree_multimap.h"
#include <iostream>
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

namespace {
/*
 * Basic idea:
 * The CondensingConverter aims to speed up relocate() operations by condensing the coordinate
 * space such that a typical relocate() will leave relocated entry in the same bucket.
 * This is achieved by "rounding" the coordinates such that a typical relocation distance is
 * smaller than the smallest quadrant (i.e. smaller than "1" -> remember that the PH-tree stores
 * internally only integer coordinates).
 *
 * The bucket_avg:
 * Taking this to the extreme would mean storing all entries in a single quadrant/bucket, e.g.
 * by setting all coordinates to "0". This would allow for extremely fast relocate() because all
 * entries are simply stored in one quadrant/bucket/hashset.
 * However, having too many entries in a single quadrant affect window query performance because
 * any window query would always match that one quadrant and would have to brute-force all entries.
 *
 * To avoid this we need to find a good middle ground where quadrants are large enough that a
 * typical relocation will stay in the same quadrant, but query does not have to filter
 * unreasonably many entries.
 *
 * A side effect of this rounding is that the PH-tree needs to store the coordinates of every entry
 * twice: The tree itself uses the rounded coordinates, but we also attach the precise coordinates
 * to every entry in order to return precise window queries and other queries.
 *
 * The implementation of this converter disregards the typical relocation distance and focuses
 * on the average bucket size.
 */
template <dimension_t DIM>
struct CondensingConverter : public ConverterPointBase<DIM, double, scalar_64_t> {
    explicit CondensingConverter(
        double estimated_area_len, size_t estimated_entity_count, size_t bucket_avg = 20) {
        // Let's assume uniform distribution for this estimate
        auto avg_bucket_count =
            (double)estimated_entity_count / (double)bucket_avg;  // = quadrant count
        auto quadrants_per_edge = pow(avg_bucket_count, 1. / DIM);
        multiplier_ = quadrants_per_edge / estimated_area_len;
        divider_ = 1. / multiplier_;
        // std::cout << "d=" << divider_ << "/" << multiplier_ << std::endl;
    }

    [[nodiscard]] PhPoint<DIM> pre(const PhPointD<DIM>& point) const {
        PhPoint<DIM> out;
        for (dimension_t i = 0; i < DIM; ++i) {
            out[i] = point[i] * multiplier_;
        }
        return out;
    }

    [[nodiscard]] PhPointD<DIM> post(const PhPoint<DIM>& in) const {
        PhPointD<DIM> out;
        for (dimension_t i = 0; i < DIM; ++i) {
            out[i] = ((double)in[i]) * divider_;
        }
        return out;
    }

    [[nodiscard]] auto pre_query(const PhBoxD<DIM>& query_box) const {
        return PhBox{pre(query_box.min()), pre(query_box.max())};
    }
    double divider_;
    double multiplier_;
};

template <typename PHTREE>
class IteratorCondBase {
    friend PHTREE;
    using V = typename PHTREE::ValueType;
    using V_INT = typename PHTREE::V_INT;

  public:
    explicit IteratorCondBase() noexcept : current_value_ptr_{nullptr} {}

    V& operator*() const noexcept {
        assert(current_value_ptr_);
        return const_cast<V&>(current_value_ptr_->first);
    }

    V* operator->() const noexcept {
        assert(current_value_ptr_);
        return const_cast<V*>(&current_value_ptr_->first);
    }

    friend bool operator==(
        const IteratorCondBase<PHTREE>& left, const IteratorCondBase<PHTREE>& right) noexcept {
        return left.current_value_ptr_ == right.current_value_ptr_;
    }

    friend bool operator!=(
        const IteratorCondBase<PHTREE>& left, const IteratorCondBase<PHTREE>& right) noexcept {
        return left.current_value_ptr_ != right.current_value_ptr_;
    }

  protected:
    void SetFinished() noexcept {
        current_value_ptr_ = nullptr;
    }

    // TODO point to V?
    void SetCurrentValue(const V_INT* current_value_ptr) noexcept {
        current_value_ptr_ = current_value_ptr;
    }

  private:
    const V_INT* current_value_ptr_;
};

template <typename ITERATOR_PH, typename PHTREE, typename FILTER>
class IteratorCondNormal : public IteratorCondBase<PHTREE> {
    friend PHTREE;

  public:
    template <typename FILT>
    explicit IteratorCondNormal(std::in_place_t, FILT&& filter) noexcept
    : IteratorCondBase<PHTREE>(), iter_ph_{}, filter_{std::forward<FILT>(filter)} {}

    template <typename ITER_PH, typename FILT>
    IteratorCondNormal(ITER_PH&& iter_ph, FILT&& filter) noexcept
    : IteratorCondBase<PHTREE>()
    , iter_ph_{std::forward<ITER_PH>(iter_ph)}
    , filter_{std::forward<FILT>(filter)} {
        FindNextElement();
    }

    IteratorCondNormal& operator++() noexcept {
        ++iter_ph_;
        FindNextElement();
        return *this;
    }

    IteratorCondNormal operator++(int) noexcept {
        IteratorCondNormal iterator(*this);
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
        while (!iter_ph_._is_finished()) {
            if (filter_(iter_ph_->second)) {
                this->SetCurrentValue(&(*iter_ph_));
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
class IteratorCondKnn : public IteratorCondNormal<ITERATOR_PH, PHTREE, FILTER> {
  public:
    template <typename ITER_PH, typename F>
    IteratorCondKnn(ITER_PH&& iter_ph, F&& filter) noexcept
    : IteratorCondNormal<ITER_PH, PHTREE, FILTER>(
          std::forward<ITER_PH>(iter_ph), std::forward<F>(filter)) {}

    [[nodiscard]] double distance() const noexcept {
        return this->GetIteratorOfPhTree().distance();
    }
};
}  // namespace

template <typename Key, typename V>
using PhEntryC = std::pair<V, Key>;

template <
    dimension_t DIM,
    typename V,
    typename CONVERTER = CondensingConverter<DIM>,
    typename BUCKET = b_plus_tree_hash_set<PhEntryC<typename CONVERTER::KeyExternal, V>>,
    bool POINT_KEYS = true,
    typename DEFAULT_QUERY_TYPE = QueryPoint>
class PhTreeMultiMapFast {
    // using CONVERTER = CondensingConverter<DIM>;
    using KeyInternal = typename CONVERTER::KeyInternal;
    using Key = typename CONVERTER::KeyExternal;
    static constexpr dimension_t DimInternal = CONVERTER::DimInternal;
    using V_INT = std::pair<V, Key>;
    using PHTREE = PhTreeMultiMapFast<DIM, V, CONVERTER, BUCKET, POINT_KEYS, DEFAULT_QUERY_TYPE>;
    using ValueType = V;
    using EndType = decltype(std::declval<PhTreeMultiMap<
                                 DIM,
                                 V_INT,
                                 CondensingConverter<DIM>,
                                 BUCKET,
                                 POINT_KEYS,
                                 DEFAULT_QUERY_TYPE>>()
                                 .end());

    friend PhTreeDebugHelper;
    friend IteratorCondBase<PHTREE>;

  public:
    using QueryBox = typename CONVERTER::QueryBoxExternal;

    // TODO remove defaults
    explicit PhTreeMultiMapFast(
        double estimated_area_len = 1000,
        size_t estimated_entity_count = 10000,
        size_t bucket_avg = 20)
    : tree_{CondensingConverter<DIM>(estimated_area_len, estimated_entity_count, bucket_avg)} {}

    explicit PhTreeMultiMapFast(CONVERTER converter) : tree_{converter} {}

    PhTreeMultiMapFast(const PhTreeMultiMapFast& other) = delete;
    PhTreeMultiMapFast& operator=(const PhTreeMultiMapFast& other) = delete;
    PhTreeMultiMapFast(PhTreeMultiMapFast&& other) noexcept = default;
    PhTreeMultiMapFast& operator=(PhTreeMultiMapFast&& other) noexcept = default;
    ~PhTreeMultiMapFast() noexcept = default;

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
    std::pair<V&, bool> emplace(const Key& key, Args&&... args) {
        auto x = tree_.try_emplace(key, wrap(key, std::forward<Args>(args)...));
        return {const_cast<V&>(unwrap(x.first)), x.second};
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
    std::pair<V&, bool> emplace_hint(const ITERATOR& iterator, const Key& key, Args&&... args) {
        auto x = tree_.emplace_hint(
            iterator.GetIteratorOfPhTree(), key, wrap(key, std::forward<Args>(args)...));
        return {const_cast<V&>(unwrap(x.first)), x.second};
    }

    /*
     * See std::unordered_multimap::insert().
     *
     * @return a pair consisting of the inserted value (or to the value that prevented the
     * insertion if the key/value already existed) and a bool denoting whether the insertion
     * took place.
     */
    std::pair<V&, bool> insert(const Key& key, const V& value) {
        return emplace(key, value);
    }

    /*
     * See emplace().
     */
    template <typename... Args>
    std::pair<V&, bool> try_emplace(const Key& key, Args&&... args) {
        return emplace(key, std::forward<Args>(args)...);
    }

    /*
     * See emplace_hint().
     */
    template <typename ITERATOR, typename... Args>
    std::pair<V&, bool> try_emplace(const ITERATOR& iterator, const Key& key, Args&&... args) {
        return emplace_hint(iterator, key, std::forward<Args>(args)...);
    }

    /*
     * @return '1', if a value is associated with the provided key, otherwise '0'.
     */
    size_t count(const Key& key) const {
        size_t n = 0;
        auto iter = tree_.find(key);
        while (iter != tree_.end()) {
            n += (iter->second == key);
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
    }

    /*
     * See std::unordered_multimap::find().
     *
     * @param key the key to look up
     * @return an iterator that points either to the first value associated with the key or
     * to {@code end()} if no value was found
     */
    auto find(const Key& key) const {
        return CreateIterator(tree_.find(key), [&key](const Key& k) { return k == key; });
    }

    /*
     * See std::unordered_multimap::find().
     *
     * @param key the key to look up
     * @param value the value to look up
     * @return an iterator that points either to the associated value of the key/value pair
     * or to {@code end()} if the key/value pair was found
     */
    auto find(const Key& key, const V& value) const {
        return CreateIterator(
            tree_.find(key, wrap(key, value)), [&key](const Key& k) { return k == key; });
    }

    /*
     * See std::unordered_multimap::erase(). Removes the provided key/value pair if it exists.
     *
     * @return '1' if the key/value pair was found, otherwise '0'.
     */
    size_t erase(const Key& key, const V& value) {
        return tree_.erase(key, wrap(key, value));
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
     * @param ignore_equal_keys Setting this to 'true' means that the function returns '0' if
     * old_key==new_key, regardless of whether old_key or new_key exist. This can be useful if the
     * result is not checked anyway or if there are a lot of "updates" with identical keys. If set
     * to 'false' (=default) it returns '0' only if old_key does not exist or if new_key has a
     * conflicting entry, otherwise it returns '1'.

     * @return '1' if a value was found and reinserted, otherwise '0'.
     */
    template <typename T2>
    size_t relocate(
        const Key& old_key, const Key& new_key, T2&& value, bool ignore_equal_keys = false) {
        auto& internal_tree = tree_._GetInternalTree();
        auto pair = internal_tree._find_or_create_two_mm(
            tree_.converter().pre(old_key), tree_.converter().pre(new_key), ignore_equal_keys);
        auto& iter_old = pair.first;
        auto& iter_new = pair.second;

        if (iter_old.IsEnd()) {
            return 0;
        }
        auto iter_old_value = iter_old->find(wrap(old_key, value));
        if (iter_old_value == iter_old->end()) {
            if (iter_new->empty()) {
                internal_tree.erase(iter_new);
            }
            return 0;
        }

        // Are we inserting in same node and same quadrant? Or are the keys equal?
        if (iter_old == iter_new) {
            assert(old_key == new_key);
            return 1;
        }

        assert(iter_old_value != iter_old->end());
        auto result = iter_new->emplace(std::move(*iter_old_value));
        if (!result.second) {
            return 0;
        }
        // TODO const_cast?
        const_cast<Key&>(result.first->second) = new_key;  // Update key in k/v pair

        iter_old->erase(iter_old_value);
        if (iter_old->empty()) {
            [[maybe_unused]] auto found = internal_tree.erase(iter_old);
            assert(found);
        }
        return 1;
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
     * @param ignore_equal_keys Setting this to 'true' means that the function returns '0' if
     * old_key==new_key, regardless of whether old_key or new_key exist. This can be useful if the
     * result is not checked anyway or if there are a lot of "updates" with identical keys. If set
     * to 'false' (=default) it returns '0' only if old_key does not exist or if new_key has a
     * conflicting entry, otherwise it returns '1'.
     *
     * @return the number of values that were relocated.
     */
    template <typename PREDICATE>
    size_t relocate_if(
        const Key& old_key,
        const Key& new_key,
        PREDICATE&& predicate,
        bool ignore_equal_keys = false) {
        auto& internal_tree = tree_._GetInternalTree();
        auto pair = internal_tree._find_or_create_two_mm(
            tree_.converter().pre(old_key), tree_.converter().pre(new_key), ignore_equal_keys);
        auto& iter_old = pair.first;
        auto& iter_new = pair.second;

        if (iter_old.IsEnd()) {
            assert(iter_new.IsEnd() || !iter_new->empty());  // Otherwise remove iter_new
            return 0;
        }

        // Are we inserting in same node and same quadrant? Or are the keys equal?
        if (iter_old == iter_new) {
            assert(old_key == new_key);
            return 1;
        }

        size_t n = 0;
        auto it = iter_old->begin();
        while (it != iter_old->end()) {
            if (predicate(unwrap(*it))) {
                auto result = iter_new->emplace(std::move(*it));
                if (result.second) {
                    it = iter_old->erase(it);
                    // TODO const_cast?
                    const_cast<Key&>(result.first->second) = new_key;
                    ++n;
                    continue;
                }
            }
            ++it;
        }

        if (iter_old->empty()) {
            [[maybe_unused]] auto found = internal_tree.erase(iter_old);
            assert(found);
        } else if (iter_new->empty()) {
            [[maybe_unused]] auto found = internal_tree.erase(iter_new);
            assert(found);
        }
        return n;
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
            WrapCallback<CALLBACK>{std::forward<CALLBACK>(callback)},
            WrapFilter<FILTER>{std::forward<FILTER>(filter)});
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
            query_box,
            WrapCallback<CALLBACK>{std::forward<CALLBACK>(callback)},
            WrapQueryFilter<QUERY_TYPE, FILTER>{
                std::forward<QUERY_TYPE>(query_type), std::forward<FILTER>(filter)},
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
        return CreateIterator(
            tree_.begin(WrapFilter<FILTER>{std::forward<FILTER>(filter)}), FilterTrue{});
        // [](const Key&) mutable { return true; });
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
        // TODO QUERY-FILTER!!!!
        return CreateIterator(
            tree_.begin_query(
                query_box,
                WrapFilter<FILTER>(std::forward<FILTER>(filter)),
                std::forward<QUERY_TYPE>(query_type)),
            [query_box](const Key& key) {
                return IsInRange(key, query_box.min(), query_box.max());
            });
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
        // TODO query filter / dist filter
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
        return IteratorCondNormal<EndType, PHTREE, FilterTrue>{std::in_place_t{}, FilterTrue{}};
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

    //    template <typename CALLBACK, typename FILTER = FilterNoOp>
    //    void for_each(CALLBACK&& callback, FILTER&& filter = FILTER()) const {
    //        tree_.for_each(
    //            NoOpCallback{},
    //            WrapCallbackFilter<CALLBACK, FILTER>{
    //                std::forward<CALLBACK>(callback), std::forward<FILTER>(filter), converter_});
    //    }
    //
    //    template <
    //        typename CALLBACK,
    //        typename FILTER = FilterNoOp,
    //        typename QUERY_TYPE = DEFAULT_QUERY_TYPE>
    //    void for_each(
    //        QueryBox query_box,
    //        CALLBACK&& callback,
    //        FILTER&& filter = FILTER(),
    //        QUERY_TYPE query_type = QUERY_TYPE()) const {
    //        tree_.template for_each<NoOpCallback, WrapCallbackFilter<CALLBACK, FILTER>>(
    //            query_type(converter_.pre_query(query_box)),
    //            {},
    //            {std::forward<CALLBACK>(callback), std::forward<FILTER>(filter), converter_});
    //    }
    //
    //    /*
    //     * Performs a rectangular window query. The parameters are the min and max keys which
    //     * contain the minimum respectively the maximum keys in every dimension.
    //     * @param query_box The query window.
    //     * @param query_type The type of query, such as QueryIntersect or QueryInclude
    //     * @param filter An optional filter function. The filter function allows filtering entries
    //     and
    //     * sub-nodes before they are returned or traversed. Any filter function must follow the
    //     * signature of the default 'FilterNoOp`.
    //     * @return Result iterator.
    //     */
    //    template <typename FILTER = FilterNoOp, typename QUERY_TYPE = DEFAULT_QUERY_TYPE>
    //    auto begin_query(
    //        const QueryBox& query_box,
    //        FILTER&& filter = FILTER(),
    //        QUERY_TYPE&& query_type = QUERY_TYPE()) const {
    //        return tree_.begin_query(query_box, std::forward<FILTER>(filter));
    //    }

    const auto& _GetInternalTree() const {
        return tree_._GetInternalTree();
    }

    void _CheckConsistencyExternal() const {
        tree_._CheckConsistencyExternal();
    }

  private:
    template <typename OUTER_ITER, typename FILTER>
    auto CreateIterator(OUTER_ITER&& outer_iter, FILTER&& filter) const {
        return IteratorCondNormal<OUTER_ITER, PHTREE, FILTER>(
            std::forward<OUTER_ITER>(outer_iter), std::forward<FILTER>(filter));
    }

    struct FilterTrue {
        template <typename KeyT>
        constexpr bool operator()(const KeyT&) const noexcept {
            return true;
        }
    };
    //    template <typename OUTER_ITER>
    //    auto CreateIteratorFind(OUTER_ITER&& outer_iter, const V& value, const Key& key) const {
    //        auto filter = [&key](const Key& key2) { return key2 == key; };
    //        return IteratorCondNormal<OUTER_ITER, PHTREE, decltype(filter)>(
    //            std::forward<OUTER_ITER>(outer_iter), std::move(filter));
    //    }
    //
    //    template <typename OUTER_ITER>
    //    auto CreateIterator(OUTER_ITER&& outer_iter, const Key& key) const {
    //        auto filter = [&key](const Key& key2) { return key2 == key; };
    //        return IteratorCondNormal<OUTER_ITER, PHTREE>(
    //            std::forward<OUTER_ITER>(outer_iter), std::move(filter));
    //    }

    template <typename OUTER_ITER>
    auto CreateIteratorKnn(OUTER_ITER&& outer_iter) const {
        // TODO
        Key key{};
        auto filter = [&key](const Key& key2) { return key2 == key; };
        return IteratorCondKnn<OUTER_ITER, PHTREE, decltype(filter)>(
            std::forward<OUTER_ITER>(outer_iter), std::move(filter));
    }

    template <typename... Args>
    auto wrap(const Key& key, Args&&... args) const {
        return std::make_pair(std::forward<Args>(args)..., key);
    }

    //    template <typename... Args>
    //    const V_INT wrap(const Key& key, Args&&... args) const {
    //        return std::make_pair(std::forward<Args>(args)..., key);
    //    }

    V& unwrap(V_INT& v) const {
        return v.first;
    }

    static const V& unwrap(const V_INT& v) {
        return v.first;
    }

    template <typename QUERY, typename FILTER>
    class WrapQueryFilter {
      public:
        template <typename Q, typename F>
        WrapQueryFilter(Q&& query, F&& filter)
        : query_{std::forward<Q>(query)}, filter_{std::forward<F>(filter)} {}

        [[nodiscard]] inline bool IsEntryValid(
            const KeyInternal& internal_key, const BUCKET& bucket) {
            return filter_.IsEntryValid(internal_key, bucket);
        }

        [[nodiscard]] inline bool IsBucketEntryValid(
            const KeyInternal& internal_key, const V_INT& value) {
            // TODO query
            return filter_.IsBucketEntryValid(internal_key, unwrap(value));
        }

        [[nodiscard]] inline bool IsNodeValid(const KeyInternal&, int) {
            return true;
        }

      private:
        QUERY query_;
        FILTER filter_;
    };

    struct NoOpCallback {
        constexpr void operator()(const Key&, const V_INT&) const noexcept {}
    };

    template <typename FILTER>
    class WrapFilter {
      public:
        // TODO merge with other constr.
        //        template <typename F>
        //        WrapFilter(F&& filter)
        //        : filter_{std::forward<F>(filter)},
        //        callback_{std::forward<NoOpCallback>(NoOpCallback{})} {}

        template <typename F>
        explicit WrapFilter(F&& filter) : filter_{std::forward<F>(filter)} {}

        WrapFilter(const WrapFilter& other) = default;
        WrapFilter& operator=(const WrapFilter& other) = default;
        WrapFilter(WrapFilter&& other) noexcept = default;
        WrapFilter& operator=(WrapFilter&& other) noexcept = default;
        ~WrapFilter() noexcept = default;

        [[nodiscard]] inline bool IsEntryValid(
            const KeyInternal& internal_key, const BUCKET& bucket) {
            return filter_.IsEntryValid(internal_key, bucket);
        }

        [[nodiscard]] inline bool IsBucketEntryValid(
            const KeyInternal& internal_key, const V_INT& value) {
            return filter_.IsBucketEntryValid(internal_key, unwrap(value));
        }

        [[nodiscard]] inline bool IsNodeValid(const KeyInternal&, int) {
            return true;
        }

      private:
        FILTER filter_;
    };

    template <typename CALLBACK = NoOpCallback>
    class WrapCallback {
      public:
        template <typename C>
        WrapCallback(C&& callback) : callback_{std::forward<C>(callback)} {}

        constexpr void operator()(const Key&, const V_INT& value) const noexcept {
            callback_(value.second, unwrap(value));
        }

      private:
        CALLBACK callback_;
    };

    template <typename DISTANCE>
    class WrapDistance {
      public:
        template <typename D>
        WrapDistance(D&& distance) : distance_{std::forward<D>(distance)} {}

        //        double operator()(const Key& p1, const Key& p2) const {
        //            // TODO unwrap!
        //            return distance_(a, b);
        //        };

      private:
        DISTANCE distance_;
    };

    PhTreeMultiMap<DIM, V_INT, CondensingConverter<DIM>, BUCKET, POINT_KEYS, DEFAULT_QUERY_TYPE>
        tree_;
};

// template <dimension_t DIM, typename T>
// auto CreateFastPhTreeMMD(
//     double estimated_area_len, size_t estimated_entity_count, size_t bucket_avg = 50) {
//     CondensingCoverter<DIM> conv(estimated_area_len, estimated_entity_count, bucket_avg);
//     return PhTreeMultiMapD<DIM, T>(conv);
// }

/**
 * A PH-Tree multi-map that uses (axis aligned) points as keys.
 * The points are defined with 64bit 'double' floating point coordinates.
 *
 * See 'PhTreeD' for details.
 */
template <
    dimension_t DIM,
    typename V,
    typename BUCKET = b_plus_tree_hash_set<PhEntryC<PhPointD<DIM>, V>>,
    typename CONVERTER = CondensingConverter<DIM>>
using PhTreeMultiMapD_C = PhTreeMultiMapFast<DIM, V, CONVERTER, BUCKET>;

// template <
//     dimension_t DIM,
//     typename T,
//     typename CONVERTER_BOX,
//     typename BUCKET = b_plus_tree_hash_set<T>>
// using PhTreeMultiMapBox_C =
//     PhTreeMultiMapFast<DIM, T, CONVERTER_BOX, BUCKET, false, QueryIntersect>;

/**
 * A PH-Tree multi-map that uses (axis aligned) boxes as keys.
 * The boxes are defined with 64bit 'double' floating point coordinates.
 *
 * See 'PhTreeD' for details.
 */
// template <dimension_t DIM, typename T, typename BUCKET = b_plus_tree_hash_set<T>>
// using PhTreeMultiMapBoxD_C = PhTreeMultiMapBox_C<DIM, T, BUCKET>;

}  // namespace improbable::phtree

namespace std {
template <typename V, typename Key>
struct std::hash<improbable::phtree::PhEntryC<Key, V>> {
    size_t operator()(const improbable::phtree::PhEntryC<Key, V>& x) const {
        return std::hash<V>{}(x.first);
    }
};

// template <class Type = void>
// struct equal_to : public binary_function<Type, Type, bool>
//{
//     bool operator()(const Type& Left, const Type& Right) const;
// };

// specialized transparent functor for operator==
template <typename Key, typename V>
struct equal_to<improbable::phtree::PhEntryC<Key, V>> {
    using X = improbable::phtree::PhEntryC<Key, V>;
//    template <class T, class U>
//    auto operator()(T&& Left, U&& Right) const
//        -> decltype(std::forward<V>(Left.first) == std::forward<V>(Right.first));

    auto operator()(const X& Left, const X& Right) const {
        return Left.first == Right.first;
    }
};

//
// template struct equal_to : binary_function
//{
//
//    // Declaration of the equal operation
//    bool operator() (const T& x,
//                    const T& y)
//        const
//    {
//        return x==y;
//    }
//
//    // Type of first parameter
//    typedef T first_argument_type;
//
//    // Type of second parameter
//    typedef T second_argument_type;
//
//    // The result is returned
//    // as bool type
//    typedef bool result_type;
//}
//
// template <typename V, typename Key>
//    bool operator() (const T& x, const T& y) const {return x==y;}
//    typedef T first_argument_type;
//    typedef T second_argument_type;
//    typedef bool result_type;
//};
};  // namespace std

#endif  // PHTREE_PHTREE_MULTIMAP_FAST_H
