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

#ifndef PHTREE_COMMON_FLAT_SPARSE_MAP_H
#define PHTREE_COMMON_FLAT_SPARSE_MAP_H

#include "bits.h"
#include <cassert>
#include <tuple>
#include <vector>

/*
 * PLEASE do not include this file directly, it is included via common.h.
 *
 * This file contains the sparse_map implementation, which is used in medium-dimensional nodes in
 * the PH-Tree.
 */
namespace improbable::phtree {

namespace {
template <typename T>
using PhSparseMapPair = std::pair<size_t, T>;

using index_t = std::int32_t;
}  // namespace

/*
 * The sparse_map is a flat map implementation that uses an array of *at* *most* SIZE=2^DIM.
 * The array contains a list sorted by key.
 *
 * It has O(log n) lookup and O(n) insertion/removal time complexity, space complexity is O(n).
 */
template <typename T>
class sparse_map {
  public:
    explicit sparse_map() : data_{} {
        data_.reserve(4);
    }

    [[nodiscard]] auto find(size_t key) {
        auto it = lower_bound(key);
        if (it != data_.end() && it->first == key) {
            return it;
        }
        return data_.end();
    }

    [[nodiscard]] auto find(size_t key) const {
        auto it = lower_bound(key);
        if (it != data_.end() && it->first == key) {
            return it;
        }
        return data_.end();
    }

    [[nodiscard]] auto lower_bound(size_t key) {
        return std::lower_bound(
            data_.begin(), data_.end(), key, [](PhSparseMapPair<T>& left, const size_t key) {
                return left.first < key;
            });
    }

    [[nodiscard]] auto lower_bound(size_t key) const {
        return std::lower_bound(
            data_.cbegin(), data_.cend(), key, [](const PhSparseMapPair<T>& left, const size_t key) {
                return left.first < key;
            });
    }

    [[nodiscard]] auto begin() {
        return data_.begin();
    }

    [[nodiscard]] auto begin() const {
        return cbegin();
    }

    [[nodiscard]] auto cbegin() const {
        return data_.cbegin();
    }

    [[nodiscard]] auto end() {
        return data_.end();
    }

    [[nodiscard]] auto end() const {
        return data_.end();
    }

    template <typename... Args>
    auto emplace(Args&&... args) {
        return try_emplace_base(std::forward<Args>(args)...);
    }

    template <typename... Args>
    auto try_emplace(size_t key, Args&&... args) {
        return try_emplace_base(key, std::forward<Args>(args)...);
    }

    void erase(size_t key) {
        auto it = lower_bound(key);
        if (it != end() && it->first == key) {
            data_.erase(it);
        }
    }

    void erase(const typename std::vector<PhSparseMapPair<T>>::iterator& iterator) {
        data_.erase(iterator);
    }

    [[nodiscard]] size_t size() const {
        return data_.size();
    }

  private:
    template <typename... Args>
    auto emplace_base(size_t key, Args&&... args) {
        auto it = lower_bound(key);
        if (it != end() && it->first == key) {
            return std::make_pair(it, false);
        } else {
            return std::make_pair(data_.emplace(it, key, std::forward<Args>(args)...), true);
        }
    }

    template <typename... Args>
    auto try_emplace_base(size_t key, Args&&... args) {
        auto it = lower_bound(key);
        if (it != end() && it->first == key) {
            return std::make_pair(it, false);
        } else {
            auto x = data_.emplace(
                it,
                std::piecewise_construct,
                std::forward_as_tuple(key),
                std::forward_as_tuple(std::forward<Args>(args)...));
            return std::make_pair(x, true);
        }
    }

    std::vector<PhSparseMapPair<T>> data_;
};

}  // namespace improbable::phtree

#endif  // PHTREE_COMMON_FLAT_SPARSE_MAP_H
