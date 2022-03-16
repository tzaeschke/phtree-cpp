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
 * connected nodes.
 *
 * The individual nodes have at most M entries.
 * The tree has O(log n) lookup and O(M log n) insertion/removal time complexity,
 * space complexity is O(n).
 */

template <typename T>
class b_plus_tree_map;

namespace {

// TODO change these type.
using key_t = scalar_64_t;
using pos_t = std::uint16_t;  // rename to index_t?
// remove this type
using index_t = scalar_64_t;

template <typename T>
class b_plus_tree_node;
template <typename T>
using BptNode = b_plus_tree_node<T>;

template <typename T>
struct BptEntry {
    BptEntry(index_t key, BptNode<T>* node) : node_{node}, second{}, first{key} {};
    template <typename... Args>
    BptEntry(index_t key, Args... args)
    : node_{nullptr}, second{std::forward<Args>(args)...}, first{key} {};
    BptNode<T>* node_;
    T second;
    index_t first;
};

template <typename T>
class BstIterator;

template <typename Entry>
using DataIterator = decltype(std::vector<BptEntry<Entry>>().begin());

constexpr static size_t M = 4;

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
template <typename T>
class b_plus_tree_node {
    using EntryT = BptEntry<T>;
    using NodeT = b_plus_tree_node<T>;
    using IterT = BstIterator<T>;
    using DataIteratorT = DataIterator<T>;
    using TreeT = b_plus_tree_map<T>;

    friend IterT;
    friend b_plus_tree_map<T>;

  public:
    explicit b_plus_tree_node(bool is_leaf, NodeT* parent, NodeT* prev, NodeT* next)
    : data_{}, is_leaf_{is_leaf}, parent_{parent}, prev_node_{prev}, next_node_{next} {};

    [[nodiscard]] auto find(index_t key) {
        auto it = lower_bound(key);
        if (it == data_.end()) {
            return it;
        }
        if (is_leaf_ && it->first == key) {
            return it;
        } else if (!is_leaf_ && it->first >= key) {
            return it;
        }
        return data_.end();
    }

    [[nodiscard]] auto find(index_t key) const {
        auto it = lower_bound(key);
        if (it != data_.end() && it->first == key) {
            return it;
        }
        return data_.end();
    }

    [[nodiscard]] auto lower_bound(index_t key) {
        return std::lower_bound(
            data_.begin(), data_.end(), key, [](EntryT& left, const index_t key) {
                return left.first < key;
            });
    }

    [[nodiscard]] auto lower_bound(index_t key) const {
        return std::lower_bound(
            data_.cbegin(), data_.cend(), key, [](const EntryT& left, const index_t key) {
                return left.first < key;
            });
    }

    //    [[nodiscard]] auto begin() {
    //        return data_.begin();
    //    }
    //
    //    [[nodiscard]] auto begin() const {
    //        return cbegin();
    //    }
    //
    //    [[nodiscard]] auto cbegin() const {
    //        return data_.cbegin();
    //    }
    //
    [[nodiscard]] auto end() {
        return data_.end();
    }

    [[nodiscard]] auto end() const {
        return data_.end();
    }
    //
    //    template <typename... Args>
    //    auto emplace(Args&&... args) {
    //        return try_emplace_base(std::forward<Args>(args)...);
    //    }

    template <typename... Args>
    auto try_emplace(const DataIteratorT& it, index_t key, TreeT& tree, Args&&... args) {
        assert(is_leaf_);
        assert(data_.size() <= M);
        if (data_.size() == M) {
            // overflow
            auto split_pos = M >> 1;
            auto split_key = data_[split_pos - 1].first;
            auto max_key = data_[M - 1].first;
            if (parent_ == nullptr) {
                // special root node handling
                assert(parent_ == nullptr);
                auto* new_parent = new NodeT(false, nullptr, nullptr, nullptr);
                new_parent->data_.emplace_back(max_key, this);
                tree.root_ = new_parent;
                parent_ = new_parent;
            }

            // create new node
            auto node2 = new NodeT(true, parent_, this, next_node_);
            if (next_node_ != nullptr) {
                next_node_->prev_node_ = node2;
            }
            next_node_ = node2;

            // populate new node
            // auto node2 = parent_->UpdateKeyAndAddNode(max_key, split_key, this, tree);
            // TODO ensure we are not incrementally inserting here, see
            //  https://stackoverflow.com/questions/15004517/moving-elements-from-stdvector-to-another-one
            node2->data_.insert(
                node2->data_.end(),
                std::make_move_iterator(data_.begin() + split_pos),
                std::make_move_iterator(data_.end()));
            // TODO: ensure we are not incrementally removing!
            data_.erase(data_.begin() + split_pos, data_.end());

            // Add node to parent
            parent_->UpdateKeyAndAddNode(max_key, split_key, std::max(max_key, key), node2, tree);
            // insert entry
            if (key <= split_key) {
                // parent_->UpdateKeyAndAddNode(max_key, split_key, this, tree);
                return EmplaceNoCheck(it, key, std::forward<Args>(args)...);
            }
            // TODO Optimize populating new node: move 1st part, insert new value, move 2nd part...?
            return node2->EmplaceNoCheck(node2->lower_bound(key), key, std::forward<Args>(args)...);
        }
        if (parent_ != nullptr) {
            auto max_key = data_[data_.size() - 1].first;
            if (key > max_key) {
                parent_->UpdateKey(max_key, key);
            }
        }
        return EmplaceNoCheck(it, key, std::forward<Args>(args)...);
    }

    // return 'true' iff an entry was erased
    bool erase(index_t key) {
        assert(is_leaf_);
        auto it = lower_bound(key);
        if (it != data_.end() && it->first == key) {
            data_.erase(it);
            return true;
        }
        return false;
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

    void _check(
        size_t& count, NodeT* parent, NodeT*& prev_leaf, key_t& known_min, key_t known_max) {
        assert(parent_ == nullptr || data_.size() >= 2);
        assert(parent_ == parent);
        if (data_.empty()) {
            assert(parent == nullptr);
            return;
        }
        assert(count == 0 || known_min < data_[0].first);
        assert(parent_ == nullptr || known_max == data_[data_.size() - 1].first);
        if (is_leaf_) {
            assert(prev_leaf == prev_node_);
            auto prev_key = data_[0].first;
            int n = 0;
            for (auto& e : data_) {
                assert(n == 0 || e.first > prev_key);
                assert(parent_ == nullptr || e.first <= known_max);
                ++n;
                ++count;
                known_min = e.first;
            }
            prev_leaf = this;
        } else {
            auto prev_key = data_[0].first;
            int n = 0;
            for (auto& e : data_) {
                assert(n == 0 || e.first > prev_key);
                e.node_->_check(count, this, prev_leaf, known_min, e.first);
                assert(parent_ == nullptr || e.first <= known_max);
                ++n;
            }
        }
    }

  private:
    void UpdateKey(index_t old_key, index_t new_key) {
        assert(new_key > old_key);
        auto it = lower_bound(old_key);
        assert(it != data_.end());
        // TODO This should work !!!
        // TODO This should work !!!
        // TODO This should work !!!
        // TODO This should work !!!
        // assert(it->first == old_key);
        assert(it->first == old_key || it->first == new_key);  // TODO remove!!!
        it->first = new_key;
        if (parent_ != nullptr && ++it == data_.end()) {
            parent_->UpdateKey(old_key, new_key);
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
        index_t key1_old, index_t key1_new, key_t key2, NodeT* child2, TreeT& tree) {
        assert(key2 > key1_new);
        assert(key1_old >= key1_new);
        assert(!is_leaf());
        NodeT* dest = this;
        bool split = false;
        if (data_.size() >= M) {
            split = true;
            auto split_pos = M >> 1;
            auto split_key = data_[split_pos - 1].first;
            auto max_key = data_[M - 1].first;
            if (parent_ == nullptr) {
                auto* new_parent = new NodeT(false, nullptr, nullptr, nullptr);
                new_parent->data_.emplace_back(max_key, this);
                tree.root_ = new_parent;
                parent_ = new_parent;
            }

            // create new node
            auto sibbling2 = new NodeT(false, parent_, this, next_node_);
            // TODO skip this for non-leaf nodes...?
            if (next_node_ != nullptr) {
                next_node_->prev_node_ = sibbling2;
            }
            next_node_ = sibbling2;

            // populate new node
            // TODO ensure we are not incrementally inserting here, see
            //  https://stackoverflow.com/questions/15004517/moving-elements-from-stdvector-to-another-one
            sibbling2->data_.insert(
                sibbling2->data_.end(),
                std::make_move_iterator(data_.begin() + split_pos),
                std::make_move_iterator(data_.end()));
            // TODO: ensure we are not incrementally removing!
            data_.erase(data_.begin() + split_pos, data_.end());
            for (auto& e : sibbling2->data_) {
                e.node_->parent_ = sibbling2;
            }

            // Add node to parent
            parent_->UpdateKeyAndAddNode(
                max_key, split_key, std::max(max_key, key2), sibbling2, tree);

            // insert entry
            if (key2 > split_key) {
                dest = sibbling2;
                child2->parent_ = dest;
                if (key1_old <= split_key) {
                    assert(data_.size() < M);
                    auto it = lower_bound(key1_old);
                    assert(it != data_.end());
                    it->first = key1_new;
                    assert(++it == data_.end());
                    assert(dest->lower_bound(key2) == dest->data_.begin());
                    dest->data_.emplace(dest->data_.begin(), key2, child2);
                    return;
                }
            }
        }
        assert(dest->data_.size() < M);
        auto it = dest->lower_bound(key1_old);
        assert(it != dest->data_.end());
        it->first = key1_new;
        ++it;
        dest->data_.emplace(it, key2, child2);
        assert(child2->parent_ == dest);
        if (parent_ != nullptr && key2 > key1_old) {
            parent_->UpdateKey(key1_old, key2);
        }
    }

    template <typename... Args>
    auto emplace_base(index_t key, Args&&... args) {
        auto it = lower_bound(key);
        if (it != data_.end() && it->first == key) {
            return std::make_pair(it, false);
        } else {
            return std::make_pair(data_.emplace(it, key, std::forward<Args>(args)...), true);
        }
    }

    template <typename... Args>
    auto try_emplace_base(const DataIteratorT& it, index_t key, Args&&... args) {
        if (it != data_.end() && it->first == key) {
            return std::make_pair(it, false);
        } else {
            //            auto x = data_.emplace(
            //                it,
            //                std::piecewise_construct,
            //                std::forward_as_tuple(key),
            //                std::forward_as_tuple(std::forward<Args>(args)...));
            auto x = data_.emplace(it, key, std::forward<Args>(args)...);
            return std::make_pair(x, true);
        }
    }

    template <typename... Args>
    auto EmplaceNoCheck(const DataIteratorT& it, index_t key, Args&&... args) {
        assert(it == data_.end() || it->first != key);
        //        if (parent_ != nullptr && it == data_.end()) {
        //            parent_->UpdateKey(data_[data_.size() - 1].first, key);
        //        }
        //            auto x = data_.emplace(
        //                it,
        //                std::piecewise_construct,
        //                std::forward_as_tuple(key),
        //                std::forward_as_tuple(std::forward<Args>(args)...));
        auto x = data_.emplace(it, key, std::forward<Args>(args)...);
        return std::make_pair(x, true);
    }

    std::vector<EntryT> data_;
    const bool is_leaf_;
    NodeT* parent_;
    NodeT* prev_node_;
    NodeT* next_node_;
};
}  // namespace

template <typename T>
class b_plus_tree_map {
    using IterT = BstIterator<T>;
    using NodeT = b_plus_tree_node<T>;

    friend NodeT;

  public:
    explicit b_plus_tree_map() : root_{new NodeT(true, nullptr, nullptr, nullptr)}, size_{0} {};

    [[nodiscard]] auto find(index_t key) {
        auto node = root_;
        while (!node->is_leaf()) {
            auto it = node->find(key);
            if (it == node->end()) {
                return end();
            }
            node = it->node_;
        }
        auto it = node->lower_bound(key);
        if (it != node->end() && it->first == key) {
            // TODO incomplete iterator, build stack?
            return IterT(node, it - node->data_.begin());
        }
        return end();
    }

    [[nodiscard]] auto find(index_t key) const {
        auto node = root_;
        while (!node->is_leaf()) {
            auto it = node->find(key);
            if (it == node->end()) {
                return end();
            }
            node = it->second.node_;
        }
        auto it = node->lower_bound(key);
        if (it != node->end() && it->first == key) {
            return it;
        }
        return end();
        //        auto it = lower_bound(key);
        //        if (it != data_.end() && it->first == key) {
        //            return it;
        //        }
        //        return data_.end();
    }

    [[nodiscard]] auto lower_bound(index_t key) {
        auto node = root_;
        while (!node->is_leaf()) {
            auto it = node->lower_bound(key);
            if (it == node->end()) {
                return end();
            }
            node = it->second.node_;
        }
        auto it = node->lower_bound(key);
        if (it != node->end() && it->first == key) {
            return it;
        }
        return end();
        //        return std::lower_bound(
        //            data_.begin(), data_.end(), key, [](PhFlatMapPair<T>& left, const size_t key)
        //            {
        //                return left.first < key;
        //            });
    }

    [[nodiscard]] auto lower_bound(index_t key) const {
        auto node = root_;
        while (!node->is_leaf()) {
            auto it = node->lower_bound(key);
            if (it == node->end()) {
                return end();
            }
            node = it->second.node_;
        }
        auto it = node->lower_bound(key);
        if (it != node->end() && it->first == key) {
            return it;
        }
        return end();
        //        return std::lower_bound(
        //            data_.cbegin(), data_.cend(), key, [](const PhFlatMapPair<T>& left, const
        //            size_t key) {
        //                return left.first < key;
        //            });
    }

    [[nodiscard]] auto begin() {
        return IterT(root_);
    }

    [[nodiscard]] auto begin() const {
        return IterT(root_);
    }

    [[nodiscard]] auto cbegin() const {
        return IterT(root_);
    }

    [[nodiscard]] auto end() {
        return IterT();
    }

    [[nodiscard]] auto end() const {
        return IterT();
    }

    template <typename... Args>
    auto emplace(Args&&... args) {
        return try_emplace_base(std::forward<Args>(args)...);
    }

    template <typename... Args>
    auto try_emplace(index_t key, Args&&... args) {
        return try_emplace_base(key, std::forward<Args>(args)...);
    }

    void erase(index_t key) {
        // TODO
        auto it = find(key);

        //        auto it = lower_bound(key);
        //        if (it != end() && it->first == key) {
        //            data_.erase(it);
        //        }
    }

    void erase(const IterT& iterator) {
        iterator.node_.erase(iterator.pos_);
    }

    [[nodiscard]] size_t size() const {
        return size_;
    }

    void _check() {
        size_t count = 0;
        NodeT* prev_leaf = nullptr;
        key_t known_min = 3423412;
        root_->_check(count, nullptr, prev_leaf, known_min, 0);
        assert(count == size());
    }

  private:
    //    template <typename... Args>
    //    auto emplace_base(size_t key, Args&&... args) {
    //        auto it = lower_bound(key);
    //        if (it != end() && it->first == key) {
    //            return std::make_pair(it, false);
    //        } else {
    //            return std::make_pair(data_.emplace(it, key, std::forward<Args>(args)...), true);
    //        }
    //    }

    template <typename... Args>
    auto try_emplace_base(key_t key, Args&&... args) {
        auto node = root_;
        while (!node->is_leaf()) {
            auto it = node->lower_bound(key);
            if (it == node->end()) {
                // insert into last node
                --it;
            }
            assert(it != node->end());
            node = it->node_;
        }
        // TODO move code into try_emplace(), or use find() ?
        auto it = node->lower_bound(key);
        if (it != node->end() && it->first == key) {
            return std::make_pair(it, false);
        }
        ++size_;
        // auto x = node->try_emplace(it, key, *this, std::forward<Args>(args)...);
        // return std::make_pair(x, true);
        return node->try_emplace(it, key, *this, std::forward<Args>(args)...);
        //        auto it = lower_bound(key);
        //        if (it != end() && it->first == key) {
        //            return std::make_pair(it, false);
        //        } else {
        //            auto x = data_.emplace(
        //                it,
        //                std::piecewise_construct,
        //                std::forward_as_tuple(key),
        //                std::forward_as_tuple(std::forward<Args>(args)...));
        //            return std::make_pair(x, true);
        //        }
    }

    //    template <typename... Args>
    //    auto try_emplace_base_inner(size_t key, Args&&... args) {
    //        auto it = lower_bound(key);
    //        if (it != end() && it->first == key) {
    //            return std::make_pair(it, false);
    //        } else {
    //            auto x = data_.emplace(
    //                it,
    //                std::piecewise_construct,
    //                std::forward_as_tuple(key),
    //                std::forward_as_tuple(std::forward<Args>(args)...));
    //            return std::make_pair(x, true);
    //        }
    //    }
    //
    //    template <typename... Args>
    //    auto try_emplace_base_leaf(size_t key, Args&&... args) {
    //        auto it = lower_bound(key);
    //        if (it != end() && it->first == key) {
    //            return std::make_pair(it, false);
    //        } else {
    //            auto x = data_.emplace(
    //                it,
    //                std::piecewise_construct,
    //                std::forward_as_tuple(key),
    //                std::forward_as_tuple(std::forward<Args>(args)...));
    //            return std::make_pair(x, true);
    //        }
    //    }

    NodeT* root_;
    size_t size_;
};

namespace {
template <typename T>
class BstIterator {
    using IterT = BstIterator<T>;
    using NodeT = b_plus_tree_node<T>;
    using EntryT = BptEntry<T>;

    friend b_plus_tree_map<T>;

  public:
    // TODO do we need this?
    // Arbitrary position iterator
    explicit BstIterator(NodeT* node, pos_t pos) : node_{node}, pos_{pos} {
        assert(node->is_leaf_ && "just for consistency, insist that we iterate leaves only ");
    }

    // begin() iterator
    explicit BstIterator(NodeT* node) : node_{node}, pos_{0} {
        assert(node->parent_ == nullptr && "must start with root node");
        // move iterator to first value
        while (!node_->is_leaf_) {
            node_ = node_->data_[0].node_;
        }
    }

    // end() iterator
    explicit BstIterator() : node_{nullptr}, pos_{0} {}

    auto& operator*() const {
        assert(AssertNotEnd());
        // TODO store pointer to entry?
        return const_cast<EntryT&>(node_->data_[pos_]);
    }

    auto* operator->() const {
        assert(AssertNotEnd());
        return const_cast<EntryT*>(&node_->data_[pos_]);
    }

    auto& operator++() {
        assert(AssertNotEnd());
        ++pos_;
        if (pos_ >= node_->data_.size()) {
            pos_ = 0;
            // this may be a nullptr -> end of data
            node_ = node_->next_node_;
        }
        return *this;
    }

    auto operator++(int) {
        IterT iterator(*this);
        ++(*this);
        return iterator;
    }

    friend bool operator==(const IterT& left, const IterT& right) {
        return left.node_ == right.node_ && left.pos_ == right.pos_;
    }

    friend bool operator!=(const IterT& left, const IterT& right) {
        return !(left == right);
    }

  private:
    [[nodiscard]] bool AssertNotEnd() const {
        return node_ != nullptr;
    }
    NodeT* node_;
    pos_t pos_;
};

template <typename T>
class BstIterator_Stack {
    using IterT = BstIterator<T>;
    using NodeT = b_plus_tree_node<T>;
    using EntryT = BptEntry<T>;

    friend b_plus_tree_map<T>;

    struct BstStackEntry {
        BstStackEntry(NodeT* node, pos_t pos) : node_{node}, pos_{pos} {};

        bool operator==(const BstStackEntry& other) const {
            return node_ == other.node_ && pos_ == other.pos_;
        }

        NodeT* node_;
        pos_t pos_;
    };

  public:
    // node=nullptr indicates end()
    explicit BstIterator_Stack(NodeT* node, pos_t pos = 0) {
        if (node != nullptr) {
            Push(node, pos);
        }
    }

    auto& operator*() const {
        assert(AssertNotEnd());
        return const_cast<EntryT&>(PeekValue());
    }

    auto* operator->() const {
        assert(AssertNotEnd());
        return const_cast<EntryT*>(&PeekValue());
    }

    auto& operator++() {
        assert(AssertNotEnd());
        auto* se = &Peek();
        while (true) {
            while (!se->node_->is_leaf() && se->node_->size() < se->pos_) {
                ++se->pos_;
                se = &Push(se->node_->data_[se->pos_].node_);
            }
            if (se->node_->size() < se->pos_) {
                ++se->pos_;
                return *this;
            }
            se = Pop();
            if (se == nullptr) {
                // finished
                return *this;
            }
            ++se->pos_;
        }
    }

    auto operator++(int) {
        IterT iterator(*this);
        ++(*this);
        return iterator;
    }

    friend bool operator==(const IterT& left, const IterT& right) {
        return left.stack_size_ == right.stack_size_ &&
            (left.stack_size_ == 0 || (left.Peek() == right.Peek()));
    }

    friend bool operator!=(const IterT& left, const IterT& right) {
        return !(left == right);
    }

  private:
    const EntryT& PeekValue() const {
        assert(stack_size_ > 0);
        auto& stack_entry = stack_[stack_size_ - 1];
        return stack_entry.node_->data_[stack_entry.pos_];
    }

    const BstStackEntry& Peek() const {
        assert(stack_size_ > 0);
        return stack_[stack_size_ - 1];
    }

    BstStackEntry& Peek() {
        assert(stack_size_ > 0);
        return stack_[stack_size_ - 1];
    }

    BstStackEntry* Pop() {
        assert(stack_size_ > 0);
        --stack_size_;
        if (stack_size_ > 0) {
            return &stack_[stack_size_ - 1];
        }
        return nullptr;
    }

    BstStackEntry& Push(NodeT* node, pos_t pos = 0) {
        if (stack_size_ + 1 > stack_.size()) {
            ++stack_size_;
            return stack_.emplace_back(node, pos);
        }
        auto& e = stack_[stack_size_];
        e.node_ = node;
        e.pos_ = pos;
        ++stack_size_;
        return e;
    }

    bool AssertNotEnd() const {
        return stack_size_ > 0;
    }
    std::vector<BstStackEntry> stack_{};
    size_t stack_size_{0};
};
}  // namespace

};  // namespace improbable::phtree

#endif  // PHTREE_COMMON_B_PLUS_TREE_H
