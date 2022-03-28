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

#ifndef PHTREE_COMMON_B_PLUS_TREE_H
#define PHTREE_COMMON_B_PLUS_TREE_H

#include "bits.h"
#include <cassert>
#include <iostream>
#include <tuple>
#include <vector>

/*
 * PLEASE do not include this file directly, it is included via common.h.
 *
 * This file contains the B+tree implementation which is used in high-dimensional nodes in
 * the PH-Tree.
 */
namespace improbable::phtree {

/*
 * The b_plus_tree_map is a B+tree implementation that uses a hierarchy of horizontally
 * connected nodes for fast traversal through all entries.
 *
 * The individual nodes have at most M entries.
 * The tree has O(log n) lookup and O(M log n) insertion/removal time complexity,
 * space complexity is O(n).
 *
 * Tree structure:
 * - Every node is either a leaf (l-node; contains values) or an inner node
 *   (n-node; contains nodes).
 * - "Sibling" nodes refer to the nodes linked by prev_node_ or next_node_. Sibling nodes
 *   usually have the same parent but may also be children of theit parent's siblings.
 *
 * - Guarantee: All leaf nodes are horizontally connected
 * - Inner nodes may or may not be connected. Specifically:
 *   - New inner nodes will be assigned siblings from the same parent or the parent's sibling
 *     (if the new node is the first or last node in a parent)
 *   - There is no guarantee that inner nodes know about their potential sibling (=other inner
 *     nodes that own bordering values/child-nodes).
 *   - There is no guarantee that siblings are on the same depth of the tree.
 * - The tree is not balanced
 */

template <typename T, size_t COUNT_MAX>
class b_plus_tree_map;

namespace {

using key_t = std::uint64_t;
using pos_t = std::uint32_t;

template <typename T, size_t COUNT_MAX>
class b_plus_tree_node;
template <typename T, size_t COUNT_MAX, typename NodeSiblingT, typename EntryT>
class b_plus_tree_node2;
template <typename T, size_t COUNT_MAX>
class b_plus_tree_node_leaf;
template <typename T, size_t COUNT_MAX>
class b_plus_tree_node_node;

template <typename T, size_t COUNT_MAX>
using BptEntryNode = std::pair<key_t, b_plus_tree_node<T, COUNT_MAX>*>;

template <typename T>
using BptEntryLeaf = std::pair<key_t, T>;

template <typename T, size_t COUNT_MAX>
class BstIterator;

/*
 * Strategy:
 * All nodes are ordered.
 *
 * Inner nodes:
 *  - have other nodes as children
 *  - the key represents the key of the value
 *
 * Leaf nodes:
 *  - have values as children.
 *  - the key represents the highest key in the child nodes
 *
 */
template <typename T, size_t COUNT_MAX>
class b_plus_tree_node {
    using NodeLT = b_plus_tree_node_leaf<T, COUNT_MAX>;
    using NodeNT = b_plus_tree_node_node<T, COUNT_MAX>;

    // TODO test with 8 again
    constexpr static size_t M_leaf = 32;
    constexpr static size_t M_inner = 32;
    constexpr static size_t M_leaf_min = 2;   // std::max((size_t)2, M_leaf >> 2); // TODO
    constexpr static size_t M_inner_min = 2;  // std::max((size_t)2, M_inner >> 2);
    // There is no point in allocating more leaf space than the max amount of entries.
    constexpr static size_t M_leaf_init = std::min(size_t(4), SIZE_MAX);
    constexpr static size_t M_inner_init = 4;

  public:
    explicit b_plus_tree_node(bool is_leaf, NodeNT* parent) noexcept
    : is_leaf_{is_leaf}, parent_{parent} {}

    virtual ~b_plus_tree_node() noexcept = default;

    template <typename E>
    [[nodiscard]] auto lower_bound(key_t key, std::vector<E>& data) noexcept {
        return std::lower_bound(data.begin(), data.end(), key, [](E& left, const key_t key) {
            return left.first < key;
        });
    }

    [[nodiscard]] inline bool is_leaf() const noexcept {
        return is_leaf_;
    }

    virtual void _check(
        size_t& count, NodeNT* parent, NodeLT*& prev_leaf, key_t& known_min, key_t known_max) = 0;

    [[nodiscard]] inline NodeNT* as_inner() noexcept {
        assert(!is_leaf_);
        return static_cast<NodeNT*>(this);
    }

    [[nodiscard]] inline NodeLT* as_leaf() noexcept {
        assert(is_leaf_);
        return static_cast<NodeLT*>(this);
    }

    [[nodiscard]] inline size_t M_min() {
        return is_leaf_ ? M_leaf_min : M_inner_min;
    }

    [[nodiscard]] inline size_t M_max() {
        return is_leaf_ ? M_leaf : M_inner;
    }

    [[nodiscard]] inline size_t M_init() {
        return is_leaf_ ? M_leaf_init : M_inner_init;
    }

  public:
    const bool is_leaf_;
    NodeNT* parent_;
};

template <typename T, size_t COUNT_MAX, typename ThisT, typename EntryT>
class b_plus_tree_node2 : public b_plus_tree_node<T, COUNT_MAX> {
    using NodeLT = b_plus_tree_node_leaf<T, COUNT_MAX>;
    using NodeNT = b_plus_tree_node_node<T, COUNT_MAX>;
    using IterT = BstIterator<T, COUNT_MAX>;
    using TreeT = b_plus_tree_map<T, COUNT_MAX>;

    friend IterT;

  public:
    explicit b_plus_tree_node2(bool is_leaf, NodeNT* parent, ThisT* prev, ThisT* next) noexcept
    : b_plus_tree_node<T, COUNT_MAX>(is_leaf, parent), data_{}, prev_node_{prev}, next_node_{next} {
        data_.reserve(this->M_init());
    }

    virtual ~b_plus_tree_node2() noexcept = default;

    [[nodiscard]] auto lower_bound(key_t key) noexcept {
        return std::lower_bound(data_.begin(), data_.end(), key, [](EntryT& left, const key_t key) {
            return left.first < key;
        });
    }

    [[nodiscard]] size_t size() const noexcept {
        return data_.size();
    }

  public:
    void Erase(pos_t pos, TreeT& tree) {
        EraseEntry(pos, tree);
    }

    void EraseEntry(pos_t pos, TreeT& tree) {
        assert(pos < data_.size());
        auto& parent_ = this->parent_;
        key_t max_key_old = data_.back().first;

        auto it_to_erase = data_.begin() + pos;
        data_.erase(it_to_erase);
        if (parent_ == nullptr) {
            if constexpr (std::is_same_v<ThisT, NodeNT>) {
                if (data_.size() < 2) {
                    auto remaining_node = data_.begin()->second;
                    data_.begin()->second = nullptr;
                    remaining_node->parent_ = nullptr;
                    tree.root_ = remaining_node;
                    delete this;
                }
            }
            return;
        }

        if (data_.size() == 0) {
            // Nothing to merge, just remove node.
            // This should be rare, i.e. only happens when a rare 1-entry node has another entry
            // removed.
            RemoveFromSiblings();
            parent_->RemoveNode(max_key_old, tree);
            return;
        }

        // TODO can we merge with different parent now?
        // TODO If not, are we unnecessarily adjusting parents anywhere?
        if (data_.size() < this->M_min()) {
            // merge
            if (prev_node_ != nullptr && prev_node_->parent_ == parent_ &&
                prev_node_->data_.size() < this->M_max()) {
                RemoveFromSiblings();
                auto& prev_data = prev_node_->data_;
                if constexpr (std::is_same_v<ThisT, NodeLT>) {
                    prev_data.emplace_back(std::move(data_[0]));
                } else {
                    data_[0].second->parent_ = prev_node_;
                    prev_data.emplace_back(std::move(data_[0]));
                    data_[0].second = nullptr;
                }
                auto prev_node = prev_node_;  // create copy because (this) will be deleted
                parent_->RemoveNode(max_key_old, tree);
                if (prev_node->parent_ != nullptr) {
                    key_t old1 = (prev_data.end() - 2)->first;
                    key_t new1 = (prev_data.end() - 1)->first;
                    prev_node->parent_->UpdateKey(old1, new1);
                }
            } else if (
                next_node_ != nullptr && next_node_->parent_ == parent_ &&
                next_node_->data_.size() < this->M_max()) {
                RemoveFromSiblings();
                assert(next_node_->parent_ == parent_);
                auto& next_data = next_node_->data_;

                if constexpr (std::is_same_v<ThisT, NodeLT>) {
                    next_data.emplace(next_data.begin(), std::move(data_[0]));
                } else {
                    data_[0].second->parent_ = next_node_;
                    next_data.emplace(next_data.begin(), std::move(data_[0]));
                    data_[0].second = nullptr;
                }

                parent_->RemoveNode(max_key_old, tree);
            } else {
                // This node is to small! Well... .
                if (pos == data_.size()) {
                    parent_->UpdateKey(max_key_old, data_[pos - 1].first);
                }
            }
        } else if (pos == data_.size()) {
            parent_->UpdateKey(max_key_old, data_[pos - 1].first);
        }
    }

    void RemoveFromSiblings() {
        if (next_node_ != nullptr) {
            next_node_->prev_node_ = prev_node_;
        }
        if (prev_node_ != nullptr) {
            prev_node_->next_node_ = next_node_;
        }
    }

    ThisT* SplitNode(key_t key, TreeT& tree, key_t& split_key, pos_t& split_pos) {
        // overflow
        auto& parent_ = this->parent_;
        split_pos = this->M_max() >> 1;
        split_key = data_[split_pos - 1].first;
        auto max_key = data_.back().first;
        if (parent_ == nullptr) {
            // special root node handling
            assert(parent_ == nullptr);
            auto* new_parent = new NodeNT(nullptr, nullptr, nullptr);
            new_parent->emplace_back(max_key, this);
            tree.root_ = new_parent;
            parent_ = new_parent;
        }

        // create new node
        auto* node2 = new ThisT(parent_, static_cast<ThisT*>(this), next_node_);
        if (next_node_ != nullptr) {
            next_node_->prev_node_ = node2;
        }
        next_node_ = node2;

        // populate new node
        node2->data_.insert(
            node2->data_.end(),
            std::make_move_iterator(data_.begin() + split_pos),
            std::make_move_iterator(data_.end()));
        data_.erase(data_.begin() + split_pos, data_.end());

        if constexpr (std::is_same_v<ThisT, NodeNT>) {
            for (auto& e : node2->data_) {
                e.second->parent_ = node2;
            }
        }

        // Add node to parent
        parent_->UpdateKeyAndAddNode(max_key, split_key, std::max(max_key, key), node2, tree);

        // Return node for insertion of new value
        return key > split_key ? node2 : static_cast<ThisT*>(this);
    }

  protected:
    std::vector<EntryT> data_;
    ThisT* prev_node_;
    ThisT* next_node_;
};

template <typename T, size_t COUNT_MAX>
class b_plus_tree_node_leaf
: public b_plus_tree_node2<T, COUNT_MAX, b_plus_tree_node_leaf<T, COUNT_MAX>, BptEntryLeaf<T>> {
    using EntryLeafT = BptEntryLeaf<T>;
    using NodeLT = b_plus_tree_node_leaf<T, COUNT_MAX>;
    using NodeNT = b_plus_tree_node_node<T, COUNT_MAX>;
    using Node2T = b_plus_tree_node2<T, COUNT_MAX, NodeLT, BptEntryLeaf<T>>;
    using IterT = BstIterator<T, COUNT_MAX>;
    using LeafIteratorT = decltype(std::vector<BptEntryLeaf<T>>().begin());
    using TreeT = b_plus_tree_map<T, COUNT_MAX>;

    friend IterT;
    friend Node2T;

  public:
    explicit b_plus_tree_node_leaf(NodeNT* parent, NodeLT* prev, NodeLT* next) noexcept
    : Node2T(true, parent, prev, next), data_{Node2T::data_} {}

    ~b_plus_tree_node_leaf() noexcept = default;

    [[nodiscard]] IterT find(key_t key) noexcept {
        auto it = this->lower_bound(key);
        if (it != data_.end() && it->first == key) {
            return IterT(this, it - data_.begin());
        }
        return IterT();
    }

    [[nodiscard]] IterT lower_bound_to_it(key_t key) noexcept {
        auto it = this->lower_bound(key);
        if (it != data_.end()) {
            return IterT(this, it - data_.begin());
        }
        return IterT();
    }

    template <typename... Args>
    auto try_emplace(key_t key, TreeT& tree, size_t& entry_count, Args&&... args) {
        auto it = this->lower_bound(key);
        if (it != data_.end() && it->first == key) {
            return std::make_pair(it, false);
        }
        ++entry_count;

        assert(data_.size() <= this->M_max());
        auto& parent_ = this->parent_;
        if (data_.size() == this->M_max()) {
            key_t split_key;
            pos_t split_pos;
            NodeLT* dest = this->SplitNode(key, tree, split_key, split_pos);

            // insert entry
            auto it_insert = it;
            if (dest != this) {
                // TODO Optimize populating new node: move 1st part, insert new value, move 2nd
                // part...? The insertion pos in node2 can be calculated:
                auto old_pos = it - data_.begin();
                it_insert = dest->data_.begin() + old_pos - split_pos;
            }
            return dest->EmplaceNoCheck(it_insert, key, std::forward<Args>(args)...);
        }
        if (parent_ != nullptr) {
            auto max_key = data_.back().first;
            if (key > max_key) {
                parent_->UpdateKey(max_key, key);
            }
        }
        return EmplaceNoCheck(it, key, std::forward<Args>(args)...);
    }

    bool erase_key(key_t key, TreeT& tree) {
        auto it = this->lower_bound(key);
        if (it != data_.end() && it->first == key) {
            this->Erase(it - data_.begin(), tree);
            return true;
        }
        return false;
    }

    void erase_pos(pos_t pos, TreeT& tree) {
        this->Erase(pos, tree);
    }

    void _check(
        size_t& count, NodeNT* parent, NodeLT*& prev_leaf, key_t& known_min, key_t known_max) {
        (void)parent;
        (void)known_max;

        // assert(parent_ == nullptr || data_.size() >= M_min);
        assert(this->parent_ == parent);
        if (data_.empty()) {
            assert(parent == nullptr);
            return;
        }
        assert(this->parent_ == nullptr || known_max == data_.back().first);
        assert(prev_leaf == this->prev_node_);
        for (auto& e : data_) {
            assert(count == 0 || e.first > known_min);
            assert(this->parent_ == nullptr || e.first <= known_max);
            ++count;
            known_min = e.first;
        }
        prev_leaf = this;
    }

  private:
    template <typename... Args>
    auto EmplaceNoCheck(const LeafIteratorT& it, key_t key, Args&&... args) {
        assert(it == data_.end() || it->first != key);
        auto x = data_.emplace(
            it,
            std::piecewise_construct,
            std::forward_as_tuple(key),
            std::forward_as_tuple(std::forward<Args>(args)...));
        return std::make_pair(x, true);
    }

    std::vector<EntryLeafT>& data_;
};

template <typename T, size_t COUNT_MAX>
class b_plus_tree_node_node : public b_plus_tree_node2<
                                  T,
                                  COUNT_MAX,
                                  b_plus_tree_node_node<T, COUNT_MAX>,
                                  BptEntryNode<T, COUNT_MAX>> {
    using EntryNodeT = BptEntryNode<T, COUNT_MAX>;
    using NodeT = b_plus_tree_node<T, COUNT_MAX>;
    using NodeLT = b_plus_tree_node_leaf<T, COUNT_MAX>;
    using NodeNT = b_plus_tree_node_node<T, COUNT_MAX>;
    using Node2T = b_plus_tree_node2<T, COUNT_MAX, NodeNT, BptEntryNode<T, COUNT_MAX>>;
    using IterT = BstIterator<T, COUNT_MAX>;
    using InnerteratorT = decltype(std::vector<BptEntryNode<T, COUNT_MAX>>().begin());
    using TreeT = b_plus_tree_map<T, COUNT_MAX>;

    friend IterT;
    friend Node2T;

  public:
    explicit b_plus_tree_node_node(NodeNT* parent, NodeNT* prev, NodeNT* next) noexcept
    : Node2T(false, parent, prev, next), data_{Node2T::data_} {}

    ~b_plus_tree_node_node() noexcept {
        for (auto& e : data_) {
            if (e.second != nullptr) {
                if (e.second->is_leaf()) {
                    delete e.second->as_leaf();
                } else {
                    delete e.second->as_inner();
                }
            }
        }
    }

    [[nodiscard]] NodeT* find(key_t key) noexcept {
        auto it = lower_bound_n(key);
        return it != data_.end() ? it->second : nullptr;
    }

    [[nodiscard]] NodeT* find_or_last(key_t key) noexcept {
        auto it = lower_bound_n(key);
        return it != data_.end() ? it->second : data_.back().second;
    }

    [[nodiscard]] auto lower_bound_n(key_t key) noexcept {
        return this->lower_bound(key);
    }

    void emplace_back(key_t key, NodeT* node) {
        data_.emplace_back(key, node);
    }

    void _check(
        size_t& count, NodeNT* parent, NodeLT*& prev_leaf, key_t& known_min, key_t known_max) {
        (void)parent;
        (void)known_max;

        // assert(parent_ == nullptr || data_.size() >= M_min);
        assert(this->parent_ == parent);
        if (data_.empty()) {
            assert(parent == nullptr);
            return;
        }
        assert(count == 0 || known_min < data_[0].first);
        assert(this->parent_ == nullptr || known_max == data_.back().first);
        auto prev_key = data_[0].first;

        int n = 0;
        for (auto& e : data_) {
            assert(n == 0 || e.first > prev_key);
            e.second->_check(count, this, prev_leaf, known_min, e.first);
            assert(this->parent_ == nullptr || e.first <= known_max);
            prev_key = e.first;
            ++n;
        }
    }

    void UpdateKey(key_t old_key, key_t new_key) {
        assert(new_key != old_key);
        auto it = lower_bound_n(old_key);
        assert(it != data_.end());
        assert(it->first == old_key);
        it->first = new_key;
        if (this->parent_ != nullptr && ++it == data_.end()) {
            this->parent_->UpdateKey(old_key, new_key);
        }
    }

    /*
     * This method does two things:
     * - It changes the key of the node (node 1) at 'key1_old' to 'key1_new' (they may be the same).
     * - It inserts a new node (node 2) after 'new_key1' with value 'key2'
     * Invariants:
     * - Node1: key1_old >= key1_new; Node 1 vs 2: key2 > new_key1
     * - there is no other key between key1_new and key1_old
     */
    void UpdateKeyAndAddNode(
        key_t key1_old, key_t key1_new, key_t key2, NodeT* child2, TreeT& tree) {
        assert(key2 > key1_new);
        assert(key1_old >= key1_new);
        NodeNT* dest = this;
        auto& parent_ = this->parent_;
        if (data_.size() >= this->M_max()) {
            key_t split_key;
            pos_t split_pos;
            dest = this->SplitNode(key2, tree, split_key, split_pos);

            // insert entry
            if (dest != this) { // if (key2 > split_key) { TODO remove
                // dest = sibling2; // TODO remove
                child2->parent_ = dest;
                if (key1_old <= split_key) {
                    assert(data_.size() < this->M_max());
                    auto it = data_.end() - 1;
                    it->first = key1_new; // TODO don't we need to parent_->UpdateKey(key1_old, key1_new);
                    assert(dest->lower_bound_n(key2) == dest->data_.begin());
                    dest->EmplaceNoCheck(dest->data_.begin(), key2, child2);
                    return;
                }
            }
        }
        auto it = dest->lower_bound_n(key1_old);
        assert(it != dest->data_.end());
        it->first = key1_new;
        ++it;
        dest->EmplaceNoCheck(it, key2, child2);
        assert(child2->parent_ == dest);
        if (parent_ != nullptr && key2 > key1_old) {
            parent_->UpdateKey(key1_old, key2);
        }
    }

    /*
     * key_old1==key_new1 signifies that they can/will be ignored.
     */
    void RemoveNode(key_t key_remove, TreeT& tree) {
        auto it_to_erase = lower_bound_n(key_remove);
        delete it_to_erase->second;
        this->EraseEntry(it_to_erase - data_.begin(), tree);
    }

  private:
    auto EmplaceNoCheck(const InnerteratorT& it, key_t key, NodeT* child) {
        assert(it == data_.end() || it->first > key);
        data_.emplace(it, key, child);
    }

    std::vector<EntryNodeT>& data_;
};
}  // namespace

template <typename T, size_t COUNT_MAX>
class b_plus_tree_map {
    using IterT = BstIterator<T, COUNT_MAX>;
    using NodeT = b_plus_tree_node<T, COUNT_MAX>;
    using NodeLT = b_plus_tree_node_leaf<T, COUNT_MAX>;
    using NodeNT = b_plus_tree_node_node<T, COUNT_MAX>;

    friend b_plus_tree_node2<T, COUNT_MAX, NodeLT, BptEntryLeaf<T>>;
    friend b_plus_tree_node2<T, COUNT_MAX, NodeNT, BptEntryNode<T, COUNT_MAX>>;
    friend NodeLT;
    friend NodeNT;

  public:
    explicit b_plus_tree_map() : root_{new NodeLT(nullptr, nullptr, nullptr)}, size_{0} {};

    ~b_plus_tree_map() {
        if (root_->is_leaf()) {
            delete root_->as_leaf();
        } else {
            delete root_->as_inner();
        }
    }

    [[nodiscard]] auto find(key_t key) noexcept {
        auto node = root_;
        while (!node->is_leaf()) {
            node = node->as_inner()->find(key);
            if (node == nullptr) {
                return end();
            }
        }
        return node->as_leaf()->find(key);
    }

    [[nodiscard]] auto find(key_t key) const noexcept {
        return const_cast<b_plus_tree_map&>(*this).find(key);
    }

    [[nodiscard]] auto lower_bound(key_t key) noexcept {
        auto node = root_;
        while (!node->is_leaf()) {
            node = node->as_inner()->find(key);
            if (node == nullptr) {
                return end();
            }
        }
        return node->as_leaf()->lower_bound_to_it(key);
    }

    [[nodiscard]] auto begin() noexcept {
        return IterT(root_);
    }

    [[nodiscard]] auto begin() const noexcept {
        return IterT(root_);
    }

    [[nodiscard]] auto cbegin() const noexcept {
        return IterT(root_);
    }

    [[nodiscard]] auto end() noexcept {
        return IterT();
    }

    [[nodiscard]] auto end() const noexcept {
        return IterT();
    }

    template <typename... Args>
    auto emplace(Args&&... args) {
        return try_emplace_base(std::forward<Args>(args)...);
    }

    template <typename... Args>
    auto try_emplace(key_t key, Args&&... args) {
        return try_emplace_base(key, std::forward<Args>(args)...);
    }

    void erase(key_t key) {
        auto node = root_;
        while (!node->is_leaf()) {
            node = node->as_inner()->find(key);
            if (node == nullptr) {
                return;
            }
        }
        size_ -= node->as_leaf()->erase_key(key, *this);
    }

    void erase(const IterT& iterator) {
        assert(iterator != end());
        --size_;
        iterator.node_->erase_pos(iterator.pos_, *this);
    }

    [[nodiscard]] size_t size() const noexcept {
        return size_;
    }

    void _check() {
        size_t count = 0;
        NodeLT* prev_leaf = nullptr;
        key_t known_min = std::numeric_limits<key_t>::max();
        root_->_check(count, nullptr, prev_leaf, known_min, 0);
        assert(count == size());
    }

  private:
    template <typename... Args>
    auto try_emplace_base(key_t key, Args&&... args) {
        auto node = root_;
        while (!node->is_leaf()) {
            node = node->as_inner()->find_or_last(key);
        }
        return node->as_leaf()->try_emplace(key, *this, size_, std::forward<Args>(args)...);
    }

    NodeT* root_;
    size_t size_;
};

namespace {
template <typename T, size_t COUNT_MAX>
class BstIterator {
    using IterT = BstIterator<T, COUNT_MAX>;
    using NodeT = b_plus_tree_node<T, COUNT_MAX>;
    using NodeLT = b_plus_tree_node_leaf<T, COUNT_MAX>;
    using EntryT = BptEntryLeaf<T>;

    friend b_plus_tree_map<T, COUNT_MAX>;

  public:
    using iterator_category = std::forward_iterator_tag;
    using value_type = T;
    using difference_type = pos_t;
    using pointer = T*;
    using reference = T&;

    // Arbitrary position iterator
    explicit BstIterator(NodeLT* node, pos_t pos) noexcept : node_{node}, pos_{pos} {
        assert(node->is_leaf_ && "just for consistency, insist that we iterate leaves only ");
    }

    // begin() iterator
    explicit BstIterator(NodeT* node) noexcept {
        assert(node->parent_ == nullptr && "must start with root node");
        // move iterator to first value
        while (!node->is_leaf_) {
            node = node->as_inner()->data_[0].second;
        }
        node_ = node->as_leaf();
        pos_ = 0;

        if (node_->size() == 0) {
            node_ = nullptr;
            return;
        }
    }

    // end() iterator
    BstIterator() noexcept : node_{nullptr}, pos_{0} {}

    auto& operator*() const noexcept {
        assert(AssertNotEnd());
        return const_cast<EntryT&>(node_->data_[pos_]);
    }

    auto* operator->() const noexcept {
        assert(AssertNotEnd());
        return const_cast<EntryT*>(&node_->data_[pos_]);
    }

    auto& operator++() noexcept {
        assert(AssertNotEnd());
        ++pos_;
        if (pos_ >= node_->data_.size()) {
            pos_ = 0;
            // this may be a nullptr -> end of data
            node_ = node_->next_node_;
        }
        return *this;
    }

    auto operator++(int) noexcept {
        IterT iterator(*this);
        ++(*this);
        return iterator;
    }

    friend bool operator==(const IterT& left, const IterT& right) noexcept {
        return left.node_ == right.node_ && left.pos_ == right.pos_;
    }

    friend bool operator!=(const IterT& left, const IterT& right) noexcept {
        return !(left == right);
    }

  private:
    [[nodiscard]] inline bool AssertNotEnd() const noexcept {
        return node_ != nullptr;
    }
    NodeLT* node_;
    pos_t pos_;
};
}  // namespace

};  // namespace improbable::phtree

#endif  // PHTREE_COMMON_B_PLUS_TREE_H
