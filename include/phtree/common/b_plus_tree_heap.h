/*
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

#ifndef PHTREE_COMMON_B_PLUS_TREE_HEAP_H
#define PHTREE_COMMON_B_PLUS_TREE_HEAP_H

#include "b_plus_tree_base.h"
#include "b_plus_tree_multimap.h"
#include "bits.h"
#include <cassert>
#include <iostream>
#include <map>
#include <tuple>
#include <vector>

/*
 * PLEASE do not include this file directly, it is included via common.h.
 *
 * This file contains the B+tree multimap implementation which is used in high-dimensional nodes in
 * the PH-Tree.
 */
namespace phtree::bptree {

/*
 *  TODO update doc!
 * The b_plus_tree_multimap is a B+tree implementation that uses a hierarchy of horizontally
 * connected nodes for fast traversal through all entries.
 *
 * Behavior
 * ========
 * This is a multimap. It behaves just like std::multimap, minus some API functions.
 * The set/map is ordered by their key.  Entries with identical keys have no specific ordering
 * but the order is stable with respect to insertion/removal of other entries.
 *
 *
 * Rationale
 * =========
 * This implementations is optimized for small entry count, however it should
 * scale well with large entry counts.
 *
 *
 * Internals
 * =========
 * The individual nodes have at most M entries.
 * The tree has O(log n) lookup and O(M log n) insertion/removal time complexity,
 * space complexity is O(n).
 *
 * Tree structure:
 * - Inner nodes: have other nodes as children; their key of an entry represents the highest
 *   key of any subnode in that entry
 * - Leaf nodes: have values as children; their key represents the key of a key/value pair
 * - Every node is either a leaf (l-node; contains values) or an inner node
 *   (n-node; contains nodes).
 * - "Sibling" nodes refer to the nodes linked by prev_node_ or next_node_. Sibling nodes
 *   usually have the same parent but may also be children of their parent's siblings.
 *
 * - Guarantee: All leaf nodes are horizontally connected
 * - Inner nodes may or may not be connected. Specifically:
 *   - New inner nodes will be assigned siblings from the same parent or the parent's sibling
 *     (if the new node is the first or last node in a parent)
 *   - There is no guarantee that inner nodes know about their potential sibling (=other inner
 *     nodes that own bordering values/child-nodes).
 *   - There is no guarantee that siblings are on the same depth of the tree.
 * - The tree is not balanced
 *
 */

namespace detail {
template <typename Value>
struct DefaultGetKey1 {
    double operator()(const Value& v) const {
        return v.first;
    }
};
}  // namespace detail

// TODO clean this up -> double?
template <
    typename Value,
    typename Compare = std::less<double>,
    typename GetKey = detail::DefaultGetKey1<Value>>
class b_plus_tree_heap {
    using Key = decltype(GetKey{}(Value{}));
    struct SwapComp {
        Compare comp;
        // template<class T>
        bool operator()(Key const& x, Key const& y) const {
            return !comp(x, y);
        }
    };

  public:
    const Value& top() const {
        assert(!data_.empty());
        return data_.back().second;
    }

    const Value& top_max() const {
        assert(!data_.empty());
        return data_.begin()->second;
    }

    Value& top_max() {
        assert(!data_.empty());
        return data_.begin()->second;
    }

    template <typename... Args>
    void emplace(Args&&... args) {
        Value v{std::forward<Args>(args)...};
        Key key = GetKey{}(v);
        data_.emplace(key, std::move(v));
        //        data_.emplace(std::forward<Args>(args)...);
    }

    void pop() {
        assert(!data_.empty());
        data_.pop_back();
    }

    void pop_max() {
        assert(!data_.empty());
        data_.erase(data_.begin());
    }

    [[nodiscard]] bool empty() const noexcept {
        return data_.empty();
    }

    [[nodiscard]] size_t size() const noexcept {
        return data_.size();
    }

    // TODO Simple hack: just negate the key to a negative value?
    void _check() const {
        data_._check();
    }

  private:
    b_plus_tree_multimap<Key, Value, Compare> data_{};  // The heap array.
};

namespace detail {
template <typename Value>
struct DefaultGetKey2 {
    double operator()(const Value& v) const {
        return v.first;
    }
};
}  // namespace detail

// TODO clean this up -> double?
template <
    typename Value,
    typename Compare = std::less<double>,
    typename GetKey = detail::DefaultGetKey2<Value>>
class b_plus_tree_heap2 {
    using Key = decltype(GetKey{}(Value{}));  // TODO remove?!?!?

  public:
    const Value& top() const {
        return data_.rbegin()->second;
    }

    const Value& top_max() const {
        return data_.begin()->second;
    }

    Value& top_max() {
        return data_.begin()->second;
    }

    template <typename... Args>
    void emplace(Args&&... args) {
        Value v{std::forward<Args>(args)...};
        Key key = GetKey{}(v);
        data_.emplace(key, std::move(v));
    }

    void pop() {
        data_.erase(--data_.end());
    }

    void pop_max() {
        data_.erase(data_.begin());
    }

    // pop_max() + emplace()
    template <typename... Args>
    void replace_max(Key& key, Args&&... args) {
        pop_max();
        emplace(key, std::forward<Args>(args)...);
    }

    [[nodiscard]] bool empty() const noexcept {
        return data_.empty();
    }

    [[nodiscard]] size_t size() const noexcept {
        return data_.size();
    }

    // TODO  Simple hack: just negate the key to a negative value?
    void _check() {
        //       data_._check();
    }

  private:
    std::multimap<Key, Value, Compare> data_;
};

}  // namespace phtree::bptree

#endif  // PHTREE_COMMON_B_PLUS_TREE_HEAP_H
