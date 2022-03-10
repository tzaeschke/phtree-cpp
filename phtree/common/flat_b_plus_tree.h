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

#ifndef PHTREE_COMMON_B_PLUS_TREE_H
#define PHTREE_COMMON_B_PLUS_TREE_H

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

/*
 * The sparse_map is a flat map implementation that uses an array of *at* *most* SIZE=2^DIM.
 * The array contains a list sorted by key.
 *
 * It has O(log n) lookup and O(n) insertion/removal time complexity, space complexity is O(n).
 */

namespace {
template <typename T>
class BptNode;

template <typename T>
struct BptEntry {
    BptNode<T>* node_;
    std::optional<T> value_;
};

template <typename T>
using BptPair = std::pair<size_t, BptEntry<T>>;

using index_t = std::int64_t;

template <typename T>
class b_plus_tree_node {
    using EntryT = BptEntry<T>;

  public:
    explicit b_plus_tree_node(bool is_leaf) : data_{}, is_leaf_{is_leaf} {};

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
            data_.begin(), data_.end(), key, [](EntryT& left, const size_t key) {
                return left.first < key;
            });
    }

    [[nodiscard]] auto lower_bound(size_t key) const {
        return std::lower_bound(
            data_.cbegin(), data_.cend(), key, [](const EntryT& left, const size_t key) {
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

    void erase(const typename std::vector<EntryT>::iterator& iterator) {
        data_.erase(iterator);
    }

    [[nodiscard]] size_t size() const {
        return data_.size();
    }

    [[nodiscard]] bool is_leaf() const {
        return is_leaf_;
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

    std::vector<EntryT> data_;
    const bool is_leaf_;
};
}  // namespace

constexpr size_t M = 8;

template <typename T>
class b_plus_tree_map {
  public:
    explicit b_plus_tree_map() : root_{true} {};

    [[nodiscard]] auto find(size_t key) {
        auto node = root_;
        while (!node.is_leaf()) {
            auto it = node.find(key);
            if (it == node.end() {
                return end();
            }
            node = it->second;
        }
        auto it = node.lower_bound(key);
        if (it != node.end() && it->first == key) {
            return it;
        }
        return end();
    }

    [[nodiscard]] auto find(size_t key) const {
        auto node = root_;
        while (!node.is_leaf()) {
            auto it = node.find(key);
            if (it == node.end() {
                return end();
            }
            node = it->second;
        }
        auto it = node.lower_bound(key);
        if (it != node.end() && it->first == key) {
            return it;
        }
        return end();
//        auto it = lower_bound(key);
//        if (it != data_.end() && it->first == key) {
//            return it;
//        }
//        return data_.end();
    }

    [[nodiscard]] auto lower_bound(size_t key) {
        return std::lower_bound(
            data_.begin(), data_.end(), key, [](PhFlatMapPair<T>& left, const size_t key) {
                return left.first < key;
            });
    }

    [[nodiscard]] auto lower_bound(size_t key) const {
        return std::lower_bound(
            data_.cbegin(), data_.cend(), key, [](const PhFlatMapPair<T>& left, const size_t key) {
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

    void erase(const typename std::vector<PhFlatMapPair<T>>::iterator& iterator) {
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

    b_plus_tree_node<T> root_;
};

}  // namespace improbable::phtree

#endif  // PHTREE_COMMON_B_PLUS_TREE_H
