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

#ifndef PHTREE_COMMON_BPT_VECTOR_TREE_H
#define PHTREE_COMMON_BPT_VECTOR_TREE_H

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <cstring>
#include <vector>

namespace phtree::bptree::detail {

template <typename V, size_t SIZE>
class bpt_vector_tree_iterator {
  private:
    using V2 = std::remove_cv_t<V>;
    using leaf_t = std::vector<V2>;
    using parent_t = std::vector<leaf_t>;
    using leaf_iter_t = std::remove_cv_t<decltype(std::declval<leaf_t>().begin())>;
    using parent_iter_t = std::remove_cv_t<decltype(std::declval<parent_t>().begin())>;

    using normal_iterator = bpt_vector_tree_iterator<V, SIZE>;

  public:
    using iterator_category = std::random_access_iterator_tag;
    using value_type = V;
    using difference_type = std::ptrdiff_t;
    using pointer = value_type*;
    using reference = value_type&;

    bpt_vector_tree_iterator() noexcept
    : parent_{nullptr}, parent_iter_{nullptr}, leaf_iter_{nullptr} {}
    explicit bpt_vector_tree_iterator(
        parent_t* parent, parent_iter_t parent_iter, leaf_iter_t leaf_iter) noexcept
    : parent_{parent}, parent_iter_{parent_iter}, leaf_iter_{leaf_iter} {}

    reference operator*() const noexcept {
        return (*leaf_iter_);
    }

    pointer operator->() const noexcept {
        return &*leaf_iter_;
    }

    constexpr bpt_vector_tree_iterator& operator++() noexcept {
        ++leaf_iter_;
        if (leaf_iter_ == parent_iter_->end()) {
            if (parent_iter_ + 1 != parent_->end()) {
                ++parent_iter_;
                leaf_iter_ = parent_iter_->begin();
            }
        }
        return *this;
    }

    //    constexpr bpt_vector_tree_iterator operator++(int) noexcept {
    //        return bpt_vector_iterator(ptr_++);
    //    }

    constexpr bool operator<(const bpt_vector_tree_iterator<V, SIZE>& right) const noexcept {
        return parent_iter_ < right.parent_iter_ ||
            (parent_iter_ == right.parent_iter_ && leaf_iter_ < right.leaf_iter_);
    }

    friend bool operator==(
        const bpt_vector_tree_iterator<V, SIZE>& left,
        const bpt_vector_tree_iterator<V, SIZE>& right) noexcept {
        return left.leaf_iter_ == right.leaf_iter_;
    }

    friend bool operator!=(
        const bpt_vector_tree_iterator<V, SIZE>& left,
        const bpt_vector_tree_iterator<V, SIZE>& right) noexcept {
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
    operator bpt_vector_tree_iterator<const V, SIZE>() {
        return bpt_vector_tree_iterator<const V, SIZE>{parent_iter_, leaf_iter_};
    }

  private:
    parent_t* parent_;
    parent_iter_t parent_iter_;
    leaf_iter_t leaf_iter_;
};
}  // namespace phtree::bptree::detail

namespace phtree::bptree {

template <typename V, size_t SIZE = 32>
class vector_tree {
    // using node_t = std::array<T, SIZE>;
    using node_t = std::vector<V>;

  public:
  public:
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
    using iterator = detail::bpt_vector_tree_iterator<value_type, SIZE>;  // LegacyRandomAccessIterator
    using const_iterator = detail::bpt_vector_tree_iterator<const value_type, SIZE>;
    // using iterator = value_type*;  // LegacyRandomAccessIterator
    // using const_iterator = const value_type*;

  public:
    vector_tree() noexcept : data_{}, size_{0} {}

    vector_tree(const vector_tree& rhs) noexcept : data_(rhs.data_), size_{rhs.size_} {}

    vector_tree(vector_tree&& rhs) noexcept : data_{std::move(rhs.data_)}, size_{rhs.size_} {}

    vector_tree& operator=(const vector_tree& rhs) noexcept {
        data_ = rhs.data_;
        size_ = rhs.size_;
        return *this;
    }

    vector_tree& operator=(vector_tree&& rhs) noexcept {
        data_ = std::move(rhs.data_);
        size_ = rhs.size_;
        return *this;
    }

    ~vector_tree() noexcept = default;

    V& operator[](size_t index) noexcept {
        assert(index < size_);
        return node(index)[index % SIZE];
    }

    const V& operator[](size_t index) const noexcept {
        assert(index < size_);
        return node(index)[index % SIZE];
    }

    template <typename... Args>
    void emplace_back(Args&&... args) {
        reserve(++size_);
        back_node().emplace_back(std::forward<Args>(args)...);
    }

    void emplace_back(V&& x) {
        reserve(++size_);
        back_node().emplace_back(std::move(x));
    }

    void emplace_back(const V& x) {
        reserve(++size_);
        back_node().emplace_back(x);
    }

    void erase_back() {
        assert(!empty());
        --size_;
        back_node().erase(back_node().end() - 1);
        if (back_node().empty() && data_.size() > 1) {
            data_.erase(data_.end() - 1);
        }
    }

    [[nodiscard]] bool empty() const noexcept {
        return size_ == 0;
    }

    [[nodiscard]] size_t size() const noexcept {
        return size_;
    }

    void reserve(size_t index) noexcept {
        while (index > data_.size() * SIZE) {
            data_.emplace_back();
        }
    }

    size_t capacity() {
        return data_.size() * SIZE;
    }

    constexpr reference front() noexcept {
        return data_.front().front();
    }

    constexpr const_reference front() const noexcept {
        return data_.front().front();
    }

    constexpr reference back() noexcept {
        return data_.back().back();
    }

    constexpr const_reference back() const noexcept {
        return data_.back().back();
    }

    constexpr iterator begin() noexcept {
        return iterator(&data_, data_.begin(), data_.front().begin());
    }

    constexpr const_iterator begin() const noexcept {
        return const_iterator(&data_, data_.begin(), data_.front().begin());
    }

    constexpr iterator end() noexcept {
        // TODO end of empty iterator???
        return iterator(&data_, data_.end() - 1, data_.back().end());
    }

    constexpr const_iterator end() const noexcept {
        // TODO end of empty iterator???
        return const_iterator(&data_, data_.end() - 1, data_.back().end());
    }

  private:
    node_t& node(size_t index) noexcept {
        return data_[index / SIZE];
    }

    const node_t& node(size_t index) const noexcept {
        return data_[index / SIZE];
    }

    node_t& back_node() noexcept {
        return data_.back();
    }

    std::vector<node_t> data_;
    size_t size_;
};

}  // namespace phtree::bptree

#endif  // PHTREE_COMMON_BPT_VECTOR_TREE_H
