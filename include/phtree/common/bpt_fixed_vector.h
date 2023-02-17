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

#ifndef PHTREE_COMMON_BPT_FIXED_VECTOR_H
#define PHTREE_COMMON_BPT_FIXED_VECTOR_H

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <cstring>
#include <tuple>

namespace phtree::bptree::detail {

template <typename V>
class bpt_vector_iterator {
  private:
    using normal_iterator = bpt_vector_iterator<V>;

  public:
    using iterator_category = std::random_access_iterator_tag;
    using value_type = V;
    using difference_type = std::ptrdiff_t;
    using pointer = value_type*;
    using reference = value_type&;

    bpt_vector_iterator() noexcept : ptr_{nullptr} {}
    explicit bpt_vector_iterator(value_type* ptr) noexcept : ptr_{ptr} {}

    reference operator*() const noexcept {
        return *ptr_;
    }

    pointer operator->() const noexcept {
        return ptr_;
    }

    constexpr bpt_vector_iterator& operator++() noexcept {
        ++ptr_;
        return *this;
    }

    constexpr bpt_vector_iterator operator++(int) noexcept {
        return bpt_vector_iterator(ptr_++);
    }

    constexpr bool operator<(const bpt_vector_iterator<V>& right) const noexcept {
        return ptr_ < right.ptr_;
    }

    friend bool operator==(
        const bpt_vector_iterator<V>& left, const bpt_vector_iterator<V>& right) noexcept {
        return left.ptr_ == right.ptr_;
    }

    friend bool operator!=(
        const bpt_vector_iterator<V>& left, const bpt_vector_iterator<V>& right) noexcept {
        return left.ptr_ != right.ptr_;
    }

    // Bidirectional iterator requirements
    constexpr normal_iterator& operator--() noexcept {
        --ptr_;
        return *this;
    }

    constexpr normal_iterator operator--(int) noexcept {
        return normal_iterator(ptr_--);
    }

    // Random access iterator requirements
    constexpr reference operator[](difference_type n) const noexcept {
        return ptr_[n];
    }

    constexpr normal_iterator& operator+=(difference_type n) noexcept {
        ptr_ += n;
        return *this;
    }

    constexpr normal_iterator operator+(difference_type n) const noexcept {
        return normal_iterator(ptr_ + n);
    }

    constexpr normal_iterator& operator-=(difference_type n) noexcept {
        ptr_ -= n;
        return *this;
    }

    constexpr normal_iterator operator-(difference_type n) const noexcept {
        return normal_iterator(ptr_ - n);
    }

    // Other // TODO???
    constexpr auto operator-(const normal_iterator& it) const noexcept {
        return ptr_ - it.ptr_;
    }

    //    constexpr normal_iterator operator-(V* ptr) const noexcept {
    //        return normal_iterator(ptr_ - ptr);
    //    }

    // implicit conversion to const iterator
    operator bpt_vector_iterator<const V>() {
        return bpt_vector_iterator<const V>{ptr_};
    }

  private:
    V* ptr_;
};

template <typename V, size_t SIZE = 16>
class bpt_vector {
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
    // using pointer 	Allocator::pointer 	(until C++11)
    //     std::allocator_traits<Allocator>::pointer 	(since C++11)
    // using const_pointer 	Allocator::const_pointer 	(until C++11)
    //             std::allocator_traits<Allocator>::const_pointer 	(since C++11)
    // TODO LegacyContiguousIterator ?!!?!?!?
    // using iterator = bpt_vector_iterator<value_type>;  // LegacyRandomAccessIterator
    // using const_iterator = bpt_vector_iterator<const value_type>;
    using iterator = value_type*;  // LegacyRandomAccessIterator
    using const_iterator = const value_type*;

    ~bpt_vector() noexcept {
        auto b = begin();
        for (auto it = b; it < b + size_; ++it) {
            it->~V();
        }
    }

    constexpr iterator begin() noexcept {
        return to_iter(0);
    }

    constexpr const_iterator begin() const noexcept {
        return to_iter_c(0);
    }

    constexpr iterator end() noexcept {
        return to_iter(size_);
    }

    constexpr const_iterator end() const noexcept {
        return to_iter_c(size_);
    }

    reference front() noexcept {
        return data(0);
    }

    const_reference front() const noexcept {
        return data_c(0);
    }

    reference back() noexcept {
        return data(size_ - 1);
    }

    const_reference back() const noexcept {
        return data_c(size_ - 1);
    }

    reference operator[](size_t index) noexcept {
        assert(index < size_);
        return data(index);
    }

    const_reference operator[](size_t index) const noexcept {
        assert(index < size_);
        return data_c(index);
    }

    template <typename InputIterator>
    // iterator insert(const_iterator iter, std::move_iterator<const_iterator> src_begin,
    // std::move_iterator<const_iterator> src_end) {
    iterator insert(const_iterator iter, InputIterator src_begin, InputIterator src_end) {
        auto length = src_end - src_begin;
        assert(size_ + length <= SIZE);
        auto index = to_index(iter);

        std::memmove(
            &*(begin() + index + length), &*(begin() + index), sizeof(V) * (size_ - index));
        //        for (size_t i = size_ + length - 1; i > index + length - 1; --i) {
        //            data(i) = std::move(data(i - length));
        //        }

        // TODO This is really hacky, it works only if "src" is also a bpt_vector....
        //        std::memmove(&*(begin() + index), &*(src_begin), sizeof(V) + length);
        auto src = src_begin;
        iterator dst = begin() + index;
        while (src != src_end) {
            *dst = std::move(*src);
            ++src;
            ++dst;
        }
        size_ += length;
        assert(size_ <= SIZE);
        return begin() + index;
    }

    template <typename... Args>
    iterator emplace(const_iterator iter, Args&&... args) {
        assert(size_ < SIZE);
        auto index = to_index(iter);

        auto src = begin() + index;
        memmove(&*(src + 1), &*src, sizeof(V) * (size_ - index));
        //        for (size_t i = size_; i > index; --i) {
        //            data(i) = std::move(data(i - 1));
        //        }

        new (reinterpret_cast<void*>(&data_[index * sizeof(V)])) V{std::forward<Args>(args)...};
        ++size_;
        return to_iter(index);
    }

    template <typename... Args>
    reference emplace_back(Args&&... args) {
        assert(size_ < SIZE);
        new (reinterpret_cast<void*>(&data_[size_ * sizeof(V)])) V{std::forward<Args>(args)...};
        ++size_;
        return data(size_ - 1);
    }

    iterator erase(const_iterator iterator) noexcept {
        auto index = to_index(iterator);
        data(index).~V();

        auto dst = begin() + index;
        memmove(&*dst, &*(dst + 1), sizeof(V) * (size_ - index - 1));
        --size_;
        return dst;
    }

    iterator erase(const_iterator first, const_iterator last) noexcept {
        auto index_0 = to_index(first);
        auto index_n = to_index(last);
        auto length = last - first;
        auto ptr_last = &*last;
        for (auto* ptr = &*first; ptr < ptr_last; ++ptr) {
            ptr->~V();
        }
        auto dst = begin() + index_0;
        memmove(&*dst, &*dst + length, sizeof(V) * (size_ - index_n));
        //        for (size_t i = index_0; i < (size_ - length); ++i) {
        //            data(i) = std::move(data(i + length));
        //        }
        size_ -= length;
        return to_iter(index_0);  // TODO return first?
    }

    [[nodiscard]] size_t size() const noexcept {
        return size_;
    }

    [[nodiscard]] bool empty() const noexcept {
        return size_ == 0;
    }

    void reserve(size_t) noexcept {}

  private:
    size_type to_index(const_iterator& iter) const noexcept {
        return &*iter - &data_c(0);
    }

    const_iterator to_iter_c(size_t index) const noexcept {
        return iterator{&data(index)};
    }

    iterator to_iter(size_t index) noexcept {
        return iterator{&data(index)};
    }

    reference data(size_t index) noexcept {
        return *std::launder(reinterpret_cast<V*>(&data_[index * sizeof(V)]));
    }

    const_reference data_c(size_t index) const noexcept {
        return *std::launder(reinterpret_cast<const V*>(&data_[index * sizeof(V)]));
    }

    // We use an untyped array to avoid implicit calls to constructors and destructors of entries.
    alignas(V) std::byte data_[sizeof(V) * SIZE];
    size_t size_{0};
};
}  // namespace phtree::bptree::detail

#endif  // PHTREE_COMMON_BPT_FIXED_VECTOR_H
