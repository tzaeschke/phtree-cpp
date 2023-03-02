/*
 * Copyright 2023 Tilmann Zäschke
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

#ifndef PHTREE_COMMON_BPT_PRIORITY_QUEUE_H
#define PHTREE_COMMON_BPT_PRIORITY_QUEUE_H

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <cstring>
#include <deque>
#include <vector>

namespace phtree::bptree::detail {

// TODO also implement bounded_queue based on std::array...?
//   -> auto-drops highest element on overflow

/**
 *
 * @tparam V
 * @tparam SIZE
 */
 // TODO do we need Compare?
template <typename V, typename Compare>
class bpt_priority_queue_iterator {
  private:
    using V2 = std::remove_cv_t<V>;
    using leaf_t = std::vector<V2>;
    using parent_t = std::vector<leaf_t>;
    using leaf_iter_t = std::remove_cv_t<decltype(std::declval<leaf_t>().begin())>;
    using parent_iter_t = std::remove_cv_t<decltype(std::declval<parent_t>().begin())>;

    using normal_iterator = bpt_priority_queue_iterator<V, Compare>;

  public:
    using iterator_category = std::random_access_iterator_tag;
    using value_type = V;
    using difference_type = std::ptrdiff_t;
    using pointer = value_type*;
    using reference = value_type&;

    bpt_priority_queue_iterator() noexcept
    : parent_{nullptr}, parent_iter_{nullptr}, leaf_iter_{nullptr} {}
    explicit bpt_priority_queue_iterator(
        parent_t* parent, parent_iter_t parent_iter, leaf_iter_t leaf_iter) noexcept
    : parent_{parent}, parent_iter_{parent_iter}, leaf_iter_{leaf_iter} {}

    reference operator*() const noexcept {
        return (*leaf_iter_);
    }

    pointer operator->() const noexcept {
        return &*leaf_iter_;
    }

    constexpr bpt_priority_queue_iterator& operator++() noexcept {
        ++leaf_iter_;
        if (leaf_iter_ == parent_iter_->end()) {
            if (parent_iter_ + 1 != parent_->end()) {
                ++parent_iter_;
                leaf_iter_ = parent_iter_->begin();
            }
        }
        return *this;
    }

    //    constexpr bpt_priority_queue_iterator operator++(int) noexcept {
    //        return bpt_vector_iterator(ptr_++);
    //    }

    constexpr bool operator<(const bpt_priority_queue_iterator<V, Compare>& right) const noexcept {
        return parent_iter_ < right.parent_iter_ ||
            (parent_iter_ == right.parent_iter_ && leaf_iter_ < right.leaf_iter_);
    }

    friend bool operator==(
        const bpt_priority_queue_iterator<V, Compare>& left,
        const bpt_priority_queue_iterator<V, Compare>& right) noexcept {
        return left.leaf_iter_ == right.leaf_iter_;
    }

    friend bool operator!=(
        const bpt_priority_queue_iterator<V, Compare>& left,
        const bpt_priority_queue_iterator<V, Compare>& right) noexcept {
        return left.leaf_iter_ != right.leaf_iter_;
    }

    // Bidirectional iterator requirements
    constexpr normal_iterator& operator--() noexcept {
        if (leaf_iter_ == parent_iter_->begin()) {
            --parent_iter_;
            leaf_iter_ = parent_iter_->end();
        }
        --leaf_iter_;
        return *this;
    }

    //    constexpr normal_iterator operator--(int) noexcept {
    //        return normal_iterator(ptr_--);
    //    }

    // Random access iterator requirements
    //    constexpr reference operator[](difference_type n) const noexcept {
    //        return ptr_[n];
    //    }
    //
    //    constexpr normal_iterator& operator+=(difference_type n) noexcept {
    //        ptr_ += n;
    //        return *this;
    //    }
    //
    //    constexpr normal_iterator operator+(difference_type n) const noexcept {
    //        return normal_iterator(ptr_ + n);
    //    }
    //
    //    constexpr normal_iterator& operator-=(difference_type n) noexcept {
    //        ptr_ -= n;
    //        return *this;
    //    }
    //
    //    constexpr normal_iterator operator-(difference_type n) const noexcept {
    //        return normal_iterator(ptr_ - n);
    //    }
    //
    //    // Other // TODO???
    //    constexpr auto operator-(const normal_iterator& it) const noexcept {
    //        return ptr_ - it.ptr_;
    //    }

    //    constexpr normal_iterator operator-(V* ptr) const noexcept {
    //        return normal_iterator(ptr_ - ptr);
    //    }

    // implicit conversion to const iterator
    operator bpt_priority_queue_iterator<const V, Compare>() {
        return bpt_priority_queue_iterator<const V, Compare>{parent_iter_, leaf_iter_};
    }

  private:
    parent_t* parent_;
    parent_iter_t parent_iter_;
    leaf_iter_t leaf_iter_;
};

/**
 * A vector_tree acts consists of a vector of vectors.
 * The idea is that it has almost the same execution speed as vector for the following operations:
 * - Access via [] operator
 * - emplace_back
 * - erase last entry (via erase_back)
 *
 * At the same time it scales much better when the vector grows because only the parent vector needs
 * resizing.
 */
template <typename V, typename Compare = std::less<>>
class priority_queue {
    // TODO implement "Member types":  https://en.cppreference.com/w/cpp/container/vector
  public:
    // Member types
    using value_type = V;
    // using allocator_type = Allocator
    using size_type = std::size_t;
    using difference_type = std::ptrdiff_t;
    using reference = value_type&;
    using const_reference = const value_type&;
    using pointer = value_type*;  //	Allocator::pointer 	(until C++11)
    //     std::allocator_traits<Allocator>::pointer 	(since C++11)
    using const_pointer = const value_type*;  //	Allocator::const_pointer 	(until C++11)
    //             std::allocator_traits<Allocator>::const_pointer 	(since C++11)
    // TODO LegacyContiguousIterator ?!!?!?!?
    using iterator =
        detail::bpt_priority_queue_iterator<value_type, Compare>;  // LegacyRandomAccessIterator
    using const_iterator = detail::bpt_priority_queue_iterator<const value_type, Compare>;
    // using iterator = value_type*;  // LegacyRandomAccessIterator
    // using const_iterator = const value_type*;

  public:
    priority_queue() noexcept : data_{}, comp_{} {}

    priority_queue(const priority_queue& rhs) noexcept : data_(rhs.data_), comp_{} {}

    priority_queue(priority_queue&& rhs) noexcept : data_{std::move(rhs.data_)}, comp_{} {}

    // TODO use default functions?
    // TODO ommit comp_? -> check std::priority_queue
    priority_queue& operator=(const priority_queue& rhs) noexcept {
        data_ = rhs.data_;
        comp_ = rhs.comp_;
        return *this;
    }

    priority_queue& operator=(priority_queue&& rhs) noexcept {
        data_ = std::move(rhs.data_);
        comp_ = std::move(rhs.comp_);
        return *this;
    }

    ~priority_queue() noexcept = default;

    const V& top() const {
        assert(!data_.empty());
        return data_.back();
    }

    // TODO rename bottom()
    const V& top_max() const {
        assert(!data_.empty());
        return data_.front();
    }

    V& top_max() {
        assert(!data_.empty());
        return data_.front();
    }

//    template <typename... Args>
//    void emplace(Args&&... args) {
//        Value v{std::forward<Args>(args)...};
//        Key key = GetKey{}(v);
//        data_.emplace(key, std::move(v));
//        //        data_.emplace(std::forward<Args>(args)...);
//    }

    void pop() {
        assert(!data_.empty());
        data_.pop_back();
    }

    void pop_max() {
        assert(!data_.empty());
        data_.erase(data_.begin());
        //data_.pop_front();
    }

    V& operator[](size_t index) noexcept {
        assert(index < data_.size());
        return data_[index];
    }

    const V& operator[](size_t index) const noexcept {
        assert(index < data_.size());
        return data_[index];
    }

    template <typename... Args>
    void emplace(Args&&... args) {
//        data_.emplace_back(std::forward<Args>(args)...);
//        std::push_heap(data_.begin(), data_.end(), comp_);
        V v{std::forward<Args>(args)...};
        // TODO this is bad!!! We should ask for key/value separately.... and avoid "first"
        auto pos = std::lower_bound(data_.begin(), data_.end(), v, comp_);
        //data_.emplace_back(std::forward<Args>(args)...);
        data_.emplace(pos, std::move(v));
    }

    // TODO remove from other collections
//    void emplace_back(V&& x) {
//        ensure_capacity(data_.size() + 1);
//        data_.emplace_back(std::move(x));
//    }

    void emplace_back(const V& v) {
        // TODO this is bad!!! We should ask for key/value separately.... and avoid "first"
        auto pos = std::lower_bound(data_.begin(), data_.end(), v, comp_);
        //data_.emplace_back(std::forward<Args>(args)...);
        data_.emplace(pos, std::move(v));
//        ensure_capacity(data_.size() + 1);
//        data_.emplace_back(x);
    }

//    void erase_back() {
//        assert(!empty());
//        data_.erase(data_.end() - 1);
//    }

    [[nodiscard]] bool empty() const noexcept {
        return data_.empty();
    }

    [[nodiscard]] size_t size() const noexcept {
        return data_.size();
    }

    void reserve(size_t size) noexcept {
        data_.reserve(size);
    }

    constexpr reference front() noexcept {
        return data_.front();
    }

    constexpr const_reference front() const noexcept {
        return data_.front();
    }

    constexpr reference back() noexcept {
        return data_.back();
    }

    constexpr const_reference back() const noexcept {
        return data_.back();
    }

    //    constexpr iterator begin() noexcept {
    //        return iterator(&data_, data_.begin(), data_.front().begin());
    //    }
    //
    //    constexpr const_iterator begin() const noexcept {
    //        return const_iterator(&data_, data_.begin(), data_.front().begin());
    //    }
    //
    //    constexpr iterator end() noexcept {
    //        // TODO end of empty iterator???
    //        return iterator(&data_, data_.end() - 1, data_.back().end());
    //    }
    //
    //    constexpr const_iterator end() const noexcept {
    //        // TODO end of empty iterator???
    //        return const_iterator(&data_, data_.end() - 1, data_.back().end());
    //    }

  private:
    void ensure_capacity(size_t index) noexcept {
//        if (index > data_.size() * SIZE) {
//            data_.emplace_back();
//        }
    }

    std::vector<V> data_;
    //std::deque<V> data_;
    Compare comp_;
};

}  // namespace phtree::bptree::detail

#endif  // PHTREE_COMMON_BPT_PRIORITY_QUEUE_H
