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

#ifndef PHTREE_COMMON_FLAT_ARRAY_MAP_H
#define PHTREE_COMMON_FLAT_ARRAY_MAP_H

#include "bits.h"
#include <bitset>
#include <cassert>
#include <tuple>

/*
 * PLEASE do not include this file directly, it is included via common.h.
 *
 * This file contains the array_map implementation, which is used in low-dimensional nodes in the
 * PH-Tree.
 */
namespace improbable::phtree {

template <typename T, std::size_t SIZE>
class flat_array_map;

namespace detail {

template <typename T>
using flat_map_pair = std::pair<size_t, T>;

template <typename T, std::size_t SIZE>
class flat_map_iterator {
    friend flat_array_map<T, SIZE>;

  public:
    flat_map_iterator() : first{0}, map_{nullptr} {};

    explicit flat_map_iterator(size_t index, const flat_array_map<T, SIZE>& map)
    : first{index}, map_{&map} {
        assert(index <= SIZE);
    }

    auto& operator*() const {
        assert(first < SIZE && map_->occupied(first));
        return const_cast<flat_map_pair<T>&>(map_->data(first));
    }

    auto* operator->() const {
        assert(first < SIZE && map_->occupied(first));
        return const_cast<flat_map_pair<T>*>(&map_->data(first));
    }

    auto& operator++() {
        first = (first + 1) >= SIZE ? SIZE : map_->lower_bound_index(first + 1);
        return *this;
    }

    auto operator++(int) {
        flat_map_iterator it(first, *map_);
        ++(*this);
        return it;
    }

    friend bool operator==(const flat_map_iterator& left, const flat_map_iterator& right) {
        return left.first == right.first;
    }

    friend bool operator!=(const flat_map_iterator& left, const flat_map_iterator& right) {
        return left.first != right.first;
    }

  private:
    size_t first;
    const flat_array_map<T, SIZE>* map_;
};
}  // namespace detail

/*
 * The array_map is a flat map implementation that uses an array of SIZE=2^DIM. The key is
 * effectively the position in the array.
 *
 * It has O(1) insertion/removal time complexity, but O(2^DIM) space complexity, so it is best used
 * when DIM is low and/or the map is known to have a high fill ratio.
 */
template <typename T, std::size_t SIZE>
class flat_array_map {
    using map_pair = detail::flat_map_pair<T>;
    using iterator = detail::flat_map_iterator<T, SIZE>;
    friend iterator;
    using bit_string_t = std::uint64_t;
    static constexpr bit_string_t U64_ONE = bit_string_t(1);

  public:
    [[nodiscard]] auto find(size_t index) noexcept {
        return occupied(index) ? iterator{index, *this} : end();
    }

    [[nodiscard]] auto lower_bound(size_t index) const {
        size_t index2 = lower_bound_index(index);
        if (index2 < SIZE) {
            return iterator{index2, *this};
        }
        return end();
    }

    [[nodiscard]] auto begin() const {
        size_t index = CountTrailingZeros(occupancy);
        // Assert index points to a valid position or outside the map if the map is empty
        assert((size() == 0 && index >= SIZE) || occupied(index));
        return iterator{index < SIZE ? index : SIZE, *this};
    }

    [[nodiscard]] auto cbegin() const {
        size_t index = CountTrailingZeros(occupancy);
        // Assert index points to a valid position or outside the map if the map is empty
        assert((size() == 0 && index >= SIZE) || occupied(index));
        return iterator{index < SIZE ? index : SIZE, *this};
    }

    [[nodiscard]] auto end() const {
        return iterator{SIZE, *this};
    }

    ~flat_array_map() noexcept {
        if (occupancy != 0) {
            for (size_t i = 0; i < SIZE; ++i) {
                if (occupied(i)) {
                    data(i).~pair();
                }
            }
        }
    }

    [[nodiscard]] size_t size() const {
        return std::bitset<64>(occupancy).count();
    }

    template <typename... Args>
    std::pair<map_pair*, bool> try_emplace_base(size_t index, Args&&... args) {
        if (!occupied(index)) {
            new (reinterpret_cast<void*>(&data_[index])) map_pair(
                std::piecewise_construct,
                std::forward_as_tuple(index),
                std::forward_as_tuple(std::forward<Args>(args)...));
            occupied(index, true);
            return {&data(index), true};
        }
        return {&data(index), false};
    }

    bool erase(size_t index) {
        if (occupied(index)) {
            data(index).~pair();
            occupied(index, false);
            return true;
        }
        return false;
    }

    bool erase(const iterator& iterator) {
        return erase(iterator.first);
    }

  private:
    /*
     * This returns the element at the given index, which is _not_ the n'th element (for n = index).
     */
    map_pair& data(size_t index) {
        assert(occupied(index));
        return *std::launder(reinterpret_cast<map_pair*>(&data_[index]));
    }

    const map_pair& data(size_t index) const {
        assert(occupied(index));
        return *std::launder(reinterpret_cast<const map_pair*>(&data_[index]));
    }

    [[nodiscard]] size_t lower_bound_index(size_t index) const {
        assert(index < SIZE);
        size_t num_zeros = CountTrailingZeros(occupancy >> index);
        // num_zeros may be equal to SIZE if no bits remain
        return std::min(SIZE, index + num_zeros);
    }

    void occupied(size_t index, bool flag) {
        (void)flag;
        assert(index < SIZE);
        assert(occupied(index) != flag);
        // flip the bit
        occupancy ^= (U64_ONE << index);
        assert(occupied(index) == flag);
    }

    [[nodiscard]] bool occupied(size_t index) const {
        return (occupancy >> index) & U64_ONE;
    }

    // We use an untyped array to avoid implicit calls to constructors and destructors of entries.
    std::aligned_storage_t<sizeof(map_pair), alignof(map_pair)> data_[SIZE];
    bit_string_t occupancy = 0;
};

/*
 * array_map is a wrapper around flat_array_map. It introduces one layer of indirection.
 * This is useful to decouple instantiation of a node from instantiation of it's descendants
 * (the flat_array_map directly instantiates an array of descendants).
 */
template <typename T, std::size_t SIZE>
class array_map {
    static_assert(SIZE <= 64);  // or else we need to adapt 'occupancy'
    static_assert(SIZE > 0);
    using iterator = improbable::phtree::detail::flat_map_iterator<T, SIZE>;

  public:
    array_map() {
        data_ = new flat_array_map<T, SIZE>();
    }

    array_map(const array_map& other) = delete;
    array_map& operator=(const array_map& other) = delete;

    array_map(array_map&& other) noexcept : data_{other.data_} {
        other.data_ = nullptr;
    }

    array_map& operator=(array_map&& other) noexcept {
        data_ = other.data_;
        other.data_ = nullptr;
        return *this;
    }

    ~array_map() {
        delete data_;
    }

    [[nodiscard]] auto find(size_t index) noexcept {
        return data_->find(index);
    }

    [[nodiscard]] auto find(size_t key) const noexcept {
        return const_cast<array_map&>(*this).find(key);
    }

    [[nodiscard]] auto lower_bound(size_t index) const {
        return data_->lower_bound(index);
    }

    [[nodiscard]] auto begin() const {
        return data_->begin();
    }

    [[nodiscard]] iterator cbegin() const {
        return data_->cbegin();
    }

    [[nodiscard]] auto end() const {
        return data_->end();
    }

    template <typename... Args>
    auto emplace(Args&&... args) {
        return data_->try_emplace_base(std::forward<Args>(args)...);
    }

    template <typename... Args>
    auto try_emplace(size_t index, Args&&... args) {
        return data_->try_emplace_base(index, std::forward<Args>(args)...);
    }

    bool erase(size_t index) {
        return data_->erase(index);
    }

    bool erase(const iterator& iterator) {
        return data_->erase(iterator);
    }

    [[nodiscard]] size_t size() const {
        return data_->size();
    }

  private:
    flat_array_map<T, SIZE>* data_;
};

}  // namespace improbable::phtree

#endif  // PHTREE_COMMON_FLAT_ARRAY_MAP_H
