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

#ifndef PHTREE_PHTREE_MULTIMAP_H
#define PHTREE_PHTREE_MULTIMAP_H

#include "common/common.h"
#include "v16/phtree_v16.h"
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
 * Base class for the internal PH-Tree multi-map iterators.
 *
 * This base class must be distinct from the other Iterator classes because it must be agnostic of
 * the types of the fields that hold iterators. If it knew about these types then we would need
 * to provide them for the ==/!= operators, which would then make it impossible to compare
 * the generic end() iterator with any specialized iterator.
 */
template <typename PHTREE>
class IteratorBase {
    friend PHTREE;
    using T = typename PHTREE::ValueType;

  public:
    explicit IteratorBase() noexcept : current_value_ptr_{nullptr} {}

    T& operator*() const noexcept {
        assert(current_value_ptr_);
        return const_cast<T&>(*current_value_ptr_);
    }

    T* operator->() const noexcept {
        assert(current_value_ptr_);
        return const_cast<T*>(current_value_ptr_);
    }

    friend bool operator==(
        const IteratorBase<PHTREE>& left, const IteratorBase<PHTREE>& right) noexcept {
        // Note: The following compares pointers to Entry objects (actually: their values T)
        // so it should be _fast_ and return 'true' only for identical entries.
        static_assert(std::is_pointer_v<decltype(left.current_value_ptr_)>);
        return left.current_value_ptr_ == right.current_value_ptr_;
    }

    friend bool operator!=(
        const IteratorBase<PHTREE>& left, const IteratorBase<PHTREE>& right) noexcept {
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
class IteratorNormal : public IteratorBase<PHTREE> {
    friend PHTREE;
    using BucketIterType = typename PHTREE::BucketIterType;

  public:
    explicit IteratorNormal() noexcept
    : IteratorBase<PHTREE>(), iter_ph_{}, iter_bucket_{}, filter_{} {}

    // TODO Why are we passing two iterators by reference + std::move?
    // See: https://abseil.io/tips/117
    IteratorNormal(ITERATOR_PH iter_ph, BucketIterType iter_bucket, const FILTER filter) noexcept
    : IteratorBase<PHTREE>()
    , iter_ph_{std::move(iter_ph)}
    , iter_bucket_{std::move(iter_bucket)}
    , filter_{filter} {
        FindNextElement();
    }

    IteratorNormal& operator++() noexcept {
        ++iter_bucket_;
        FindNextElement();
        return *this;
    }

    IteratorNormal operator++(int) noexcept {
        IteratorNormal iterator(*this);
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
    auto& GetIteratorOfBucket() const noexcept {
        return iter_bucket_;
    }

    auto& GetIteratorOfPhTree() const noexcept {
        return iter_ph_;
    }

  private:
    void FindNextElement() {
        while (!iter_ph_.IsEnd()) {
            while (iter_bucket_ != iter_ph_->end()) {
                // We filter only entries here, nodes are filtered elsewhere
                if (filter_.IsEntryValid(iter_ph_.GetCurrentResult()->GetKey(), *iter_bucket_)) {
                    this->SetCurrentValue(&(*iter_bucket_));
                    return;
                }
                ++iter_bucket_;
            }
            ++iter_ph_;
            if (!iter_ph_.IsEnd()) {
                iter_bucket_ = iter_ph_->begin();
            }
        }
        // finished
        this->SetFinished();
    }

    ITERATOR_PH iter_ph_;
    BucketIterType iter_bucket_;
    FILTER filter_;
};

template <typename ITERATOR_PH, typename PHTREE, typename FILTER>
class IteratorKnn : public IteratorNormal<ITERATOR_PH, PHTREE, FILTER> {
    using BucketIterType = typename PHTREE::BucketIterType;

  public:
    IteratorKnn(const ITERATOR_PH iter_ph, BucketIterType iter_bucket, const FILTER filter) noexcept
    : IteratorNormal<ITERATOR_PH, PHTREE, FILTER>(iter_ph, iter_bucket, filter) {}

    [[nodiscard]] double distance() const noexcept {
        return this->GetIteratorOfPhTree().distance();
    }
};

}  // namespace

/*
 * The PhTreeMultiMap class.
 */
template <
    dimension_t DIM,
    typename T,
    typename CONVERTER = ConverterNoOp<DIM, scalar_64_t>,
    typename BUCKET = std::unordered_set<T>,
    bool POINT_KEYS = true,
    typename DEFAULT_QUERY_TYPE = QueryPoint>
class PhTreeMultiMap {
    friend PhTreeDebugHelper;
    using KeyInternal = typename CONVERTER::KeyInternal;
    using QueryBox = typename CONVERTER::QueryBoxExternal;
    using Key = typename CONVERTER::KeyExternal;
    static constexpr dimension_t DimInternal = CONVERTER::DimInternal;
    using PHTREE = PhTreeMultiMap<DIM, T, CONVERTER, BUCKET, POINT_KEYS, DEFAULT_QUERY_TYPE>;

  public:
    using ValueType = T;
    using BucketIterType = decltype(std::declval<BUCKET>().begin());
    using EndType = decltype(std::declval<v16::PhTreeV16<DimInternal, BUCKET, CONVERTER>>().end());

    explicit PhTreeMultiMap(CONVERTER converter = CONVERTER())
    : tree_{&converter_}, converter_{converter}, size_{0} {}

    PhTreeMultiMap(const PhTreeMultiMap& other) = delete;
    PhTreeMultiMap& operator=(const PhTreeMultiMap& other) = delete;
    PhTreeMultiMap(PhTreeMultiMap&& other) noexcept = default;
    PhTreeMultiMap& operator=(PhTreeMultiMap&& other) noexcept = default;
    ~PhTreeMultiMap() noexcept = default;

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
        auto& outer_iter = tree_.emplace(converter_.pre(key)).first;
        auto bucket_iter = outer_iter.emplace(std::forward<Args>(args)...);
        size_ += bucket_iter.second ? 1 : 0;
        return {const_cast<T&>(*bucket_iter.first), bucket_iter.second};
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
     * // Move value from key1 to key2
     * auto iter = tree.find(key1);
     * auto value = iter.second(); // The value may become invalid in erase()
     * erase(iter);
     * emplace_hint(iter, key2, value);  // the iterator can still be used as hint here
     */
    template <typename ITERATOR, typename... Args>
    std::pair<T&, bool> emplace_hint(const ITERATOR& iterator, const Key& key, Args&&... args) {
        auto result_ph = tree_.emplace_hint(iterator.GetIteratorOfPhTree(), converter_.pre(key));
        auto& bucket = result_ph.first;
        if (result_ph.second) {
            // new bucket
            auto result = bucket.emplace(std::forward<Args>(args)...);
            size_ += result.second;
            return {const_cast<T&>(*result.first), result.second};
        } else {
            // existing bucket -> we can use emplace_hint with iterator
            size_t old_size = bucket.size();
            auto result =
                bucket.emplace_hint(iterator.GetIteratorOfBucket(), std::forward<Args>(args)...);
            bool success = old_size < bucket.size();
            size_ += success;
            return {const_cast<T&>(*result), success};
        }
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
     * @return '1', if a value is associated with the provided key, otherwise '0'.
     */
    size_t count(const Key& key) const {
        auto iter = tree_.find(converter_.pre(key));
        if (iter != tree_.end()) {
            return iter->size();
        }
        return 0;
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
        auto counter_lambda = [&](const Key&, const BUCKET& bucket) { n += bucket.size(); };
        tree_.for_each(query_type(converter_.pre_query(query_box)), counter_lambda);
        return n;
    }

    /*
     * See std::unordered_multimap::find().
     *
     * @param key the key to look up
     * @return an iterator that points either to the the first value associated with the key or
     * to {@code end()} if no value was found
     */
    auto find(const Key& key) const {
        auto outer_iter = tree_.find(converter_.pre(key));
        if (outer_iter == tree_.end()) {
            return CreateIterator(outer_iter, bucket_dummy_end_);
        }
        auto bucket_iter = outer_iter.second().begin();
        return CreateIterator(outer_iter, bucket_iter);
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
        auto outer_iter = tree_.find(converter_.pre(key));
        if (outer_iter == tree_.end()) {
            return CreateIterator(outer_iter, bucket_dummy_end_);
        }
        auto bucket_iter = outer_iter.second().find(value);
        return CreateIterator(outer_iter, bucket_iter);
    }

    /*
     * See std::unordered_multimap::erase(). Removes the provided key/value pair if it exists.
     *
     * @return '1' if the key/value pair was found, otherwise '0'.
     */
    size_t erase(const Key& key, const T& value) {
        auto iter_outer = tree_.find(converter_.pre(key));
        if (iter_outer != tree_.end()) {
            auto& bucket = *iter_outer;
            auto result = bucket.erase(value);
            if (bucket.empty()) {
                tree_.erase(iter_outer);
            }
            size_ -= result;
            return result;
        }
        return 0;
    }

    /*
     * See std::map::erase(). Removes any entry located at the provided iterator.
     *
     * This function uses the iterator to directly erase the entry so it is usually faster than
     * erase(key, value).
     *
     * @return '1' if a value was found, otherwise '0'.
     */
    template <typename ITERATOR>
    size_t erase(const ITERATOR& iterator) {
        static_assert(
            std::is_convertible_v<ITERATOR*, IteratorBase<PHTREE>*>,
            "erase(iterator) requires an iterator argument. For erasing by key please use "
            "erase(key, value).");
        if (iterator != end()) {
            auto& bucket = const_cast<BUCKET&>(*iterator.GetIteratorOfPhTree());
            size_t old_size = bucket.size();
            bucket.erase(iterator.GetIteratorOfBucket());
            bool success = bucket.size() < old_size;
            if (bucket.empty()) {
                success &= tree_.erase(iterator.GetIteratorOfPhTree()) > 0;
            }
            size_ -= success;
            return success;
        }
        return 0;
    }

    /*
     * This function attempts to remove the 'value' from 'old_key' and reinsert it for 'new_key'.
     *
     * The relocate will report _success_ in the following cases:
     * - the value was removed from the old position and reinserted at the new position
     * - the position and new position refer to the same bucket.
     *
     * The relocate will report_failure_ in the following cases:
     * - The value was already present in the new position
     * - The value was not present in the old position
     *
     * This method will _always_ attempt to insert the value at the new position even if the value
     * was not found at the old position.
     * This method will _not_ remove the value from the old position if it is already present at the
     * new position.
     *
     * @param old_key The old position
     * @param new_key The new position
     * @param always_erase Setting this flag to 'true' ensures that the value is removed from
     * the old position even if it is already present at the new position. This may double the
     * execution cost of this method. The default is 'false'.
     * @return '1' if a value was found and reinserted, otherwise '0'.
     */
    size_t relocate(
        const Key& old_key, const Key& new_key, const T& value, bool always_erase = false) {
        // Be smart: insert first, if the target-map already contains the entry we can avoid erase()
        auto new_key_pre = converter_.pre(new_key);
        auto& new_bucket = tree_.emplace(new_key_pre).first;
        auto new_result = new_bucket.emplace(value);
        if (!new_result.second) {
            // Entry is already in correct place -> abort
            // Return '1' if old/new refer to the same bucket, otherwise '0'
            if (converter_.pre(old_key) == new_key_pre) {
                return 1;
            }
            if (!always_erase) {
                // Abort, unless we insist on erase()
                return 0;
            }
        }

        auto old_outer_iter = tree_.find(converter_.pre(old_key));
        if (old_outer_iter == tree_.end()) {
            // No entry for old_key -> fail
            return 0;
        }

        auto old_bucket_iter = old_outer_iter->find(value);
        if (old_bucket_iter == old_outer_iter->end()) {
            return 0;
        }
        old_outer_iter->erase(old_bucket_iter);

        // clean up
        if (old_outer_iter->empty()) {
            tree_.erase(old_outer_iter);
        }
        return 1;
    }

    /*
     * Iterates over all entries in the tree. The optional filter allows filtering entries and nodes
     * (=sub-trees) before returning / traversing them. By default all entries are returned. Filter
     * functions must implement the same signature as the default 'FilterNoOp'.
     *
     * @param callback The callback function to be called for every entry that matches the filter.
     * The callback requires the following signature: callback(const PhPointD<DIM> &, const T &)
     * @param filter An optional filter function. The filter function allows filtering entries and
     * sub-nodes before they are passed to the callback or traversed. Any filter function must
     * follow the signature of the default 'FilterNoOp`.
     * The default 'FilterNoOp` filter matches all entries.
     */
    template <typename CALLBACK_FN, typename FILTER = FilterNoOp>
    void for_each(CALLBACK_FN& callback, FILTER filter = FILTER()) const {
        CallbackWrapper<CALLBACK_FN, FILTER> inner_callback{callback, filter, converter_};
        tree_.for_each(inner_callback, WrapFilter(filter));
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
        typename CALLBACK_FN,
        typename FILTER = FilterNoOp,
        typename QUERY_TYPE = DEFAULT_QUERY_TYPE>
    void for_each(
        QueryBox query_box,
        CALLBACK_FN& callback,
        const FILTER& filter = FILTER(),
        QUERY_TYPE query_type = QUERY_TYPE()) const {
        CallbackWrapper<CALLBACK_FN, FILTER> inner_callback{callback, filter, converter_};
        tree_.for_each(
            query_type(converter_.pre_query(query_box)), inner_callback, WrapFilter(filter));
    }

    /*
     * Iterates over all entries in the tree. The optional filter allows filtering entries and nodes
     * (=sub-trees) before returning / traversing them. By default all entries are returned. Filter
     * functions must implement the same signature as the default 'FilterNoOp'.
     *
     * @return an iterator over all (filtered) entries in the tree,
     */
    template <typename FILTER = FilterNoOp>
    auto begin(FILTER filter = FILTER()) const {
        auto outer_iter = tree_.begin(WrapFilter(filter));
        if (outer_iter == tree_.end()) {
            return CreateIterator(outer_iter, bucket_dummy_end_, filter);
        }
        auto bucket_iter = outer_iter.second().begin();
        assert(bucket_iter != outer_iter.second().end());
        return CreateIterator(outer_iter, bucket_iter, filter);
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
        FILTER filter = FILTER(),
        QUERY_TYPE query_type = QUERY_TYPE()) const {
        auto outer_iter =
            tree_.begin_query(query_type(converter_.pre_query(query_box)), WrapFilter(filter));
        if (outer_iter == tree_.end()) {
            return CreateIterator(outer_iter, bucket_dummy_end_, filter);
        }
        auto bucket_iter = outer_iter.second().begin();
        assert(bucket_iter != outer_iter.second().end());
        return CreateIterator(outer_iter, bucket_iter, filter);
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
        DISTANCE distance_function = DISTANCE(),
        FILTER filter = FILTER()) const {
        // We use pre() instead of pre_query() here because, strictly speaking, we want to
        // find the nearest neighbors of a (fictional) key, which may as well be a box.
        auto outer_iter = tree_.begin_knn_query(
            min_results, converter_.pre(center), distance_function, WrapFilter(filter));
        if (outer_iter == tree_.end()) {
            return CreateIteratorKnn(outer_iter, bucket_dummy_end_, filter);
        }
        auto bucket_iter = outer_iter.second().begin();
        assert(bucket_iter != outer_iter.second().end());
        return CreateIteratorKnn(outer_iter, bucket_iter, filter);
    }

    /*
     * @return An iterator representing the tree's 'end'.
     */
    auto end() const {
        return IteratorNormal<EndType, PHTREE, FilterNoOp>{};
    }

    /*
     * Remove all entries from the tree.
     */
    void clear() {
        tree_.clear();
        size_ = 0;
    }

    /*
     * @return the number of entries (key/value pairs) in the tree.
     */
    [[nodiscard]] size_t size() const {
        return size_;
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

    template <typename OUTER_ITER, typename FILTER = FilterNoOp>
    auto CreateIterator(
        OUTER_ITER outer_iter, BucketIterType bucket_iter, FILTER filter = FILTER()) const {
        return IteratorNormal<OUTER_ITER, PHTREE, FILTER>(
            std::move(outer_iter), std::move(bucket_iter), filter);
    }

    template <typename OUTER_ITER, typename FILTER = FilterNoOp>
    auto CreateIteratorKnn(
        OUTER_ITER outer_iter, BucketIterType bucket_iter, FILTER filter = FILTER()) const {
        return IteratorKnn<OUTER_ITER, PHTREE, FILTER>(
            std::move(outer_iter), std::move(bucket_iter), filter);
    }

    template <typename FILTER>
    static auto WrapFilter(FILTER filter) {
        // We always have two iterators, one that traverses the PH-Tree and one that traverses the
        // bucket. Using the FilterWrapper we create a new Filter for the PH-Tree iterator. This new
        // filter checks only if nodes are valid. It cannot check whether buckets are valid.
        // The original filter is then used when we iterate over the entries of a bucket. At this
        // point, we do not need to check IsNodeValid anymore for each entry (see `IteratorNormal`).
        struct FilterWrapper {
            [[nodiscard]] constexpr bool IsEntryValid(const KeyInternal&, const BUCKET&) const {
                // This filter is checked in the Iterator.
                return true;
            }
            [[nodiscard]] constexpr bool IsNodeValid(
                const KeyInternal& prefix, int bits_to_ignore) const {
                return filter_.IsNodeValid(prefix, bits_to_ignore);
            }
            FILTER filter_;
        };
        return FilterWrapper{filter};
    }

    template <typename CALLBACK_FN, typename FILTER>
    struct CallbackWrapper {
        /*
         * The CallbackWrapper ensures that we call the callback on each entry of the bucket.
         * The vanilla PH-Tree call it only on the bucket itself.
         */
        void operator()(const Key& key, const BUCKET& bucket) const {
            auto internal_key = converter_.pre(key);
            for (auto& entry : bucket) {
                if (filter_.IsEntryValid(internal_key, entry)) {
                    callback_(key, entry);
                }
            }
        }
        CALLBACK_FN& callback_;
        const FILTER filter_;
        const CONVERTER& converter_;
    };

    v16::PhTreeV16<DimInternal, BUCKET, CONVERTER> tree_;
    CONVERTER converter_;
    BucketIterType bucket_dummy_end_;
    size_t size_;
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
    typename CONVERTER = ConverterIEEE<DIM>,
    typename BUCKET = std::unordered_set<T>>
using PhTreeMultiMapD = PhTreeMultiMap<DIM, T, CONVERTER, BUCKET>;

template <
    dimension_t DIM,
    typename T,
    typename CONVERTER_BOX,
    typename BUCKET = std::unordered_set<T>>
using PhTreeMultiMapBox = PhTreeMultiMap<DIM, T, CONVERTER_BOX, BUCKET, false, QueryIntersect>;

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
    typename BUCKET = std::unordered_set<T>>
using PhTreeMultiMapBoxD = PhTreeMultiMapBox<DIM, T, CONVERTER_BOX, BUCKET>;

}  // namespace improbable::phtree

#endif  // PHTREE_PHTREE_MULTIMAP_H
