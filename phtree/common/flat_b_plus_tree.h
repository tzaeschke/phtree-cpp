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

constexpr static size_t M = 88;

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
        if (it != data_.end() && it->first == key) {
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
            auto split_key = data_[split_pos].first;
            auto max_key = data_[M - 1].first;
            if (parent_) {
                // TODO
                assert(false && "Implement overflow");
//                auto split_pos = M >> 1;
//                auto split_key = data_[split_pos].first;
//                auto max_key = data_[M - 1].first;
                auto node2 = parent_->UpdateKeyAndAddNode(max_key, split_key, this);
                // TODO ensure we are not incrementally inserting here, see
                //  https://stackoverflow.com/questions/15004517/moving-elements-from-stdvector-to-another-one
                node2->data_.insert(node2->data_.end(), std::make_move_iterator(data_.begin() + split_pos),
                          std::make_move_iterator(data_.end()));
                // TODO: ensure we are not incrementally removing!
                data_.erase(data_.begin() + split_pos, data_.end());

                // insert entry
                return try_emplace_base(data_.end(), key, std::forward<Args>(args)...);
            } else {
                assert(parent_ == nullptr);
                auto* new_parent = new NodeT(false, nullptr, nullptr, nullptr);
                new_parent->data_.emplace_back(max_key, this);
                tree.root_ = new_parent;
                parent_ = new_parent;


                auto node2 = parent_->UpdateKeyAndAddNode(max_key, split_key, this);
                // TODO ensure we are not incrementally inserting here, see
                //  https://stackoverflow.com/questions/15004517/moving-elements-from-stdvector-to-another-one
                node2->data_.insert(node2->data_.end(), std::make_move_iterator(data_.begin() + split_pos),
                                    std::make_move_iterator(data_.end()));
                // TODO: ensure we are not incrementally removing!
                data_.erase(data_.begin() + split_pos, data_.end());

                // insert entry
                return try_emplace_base(data_.end(), key, std::forward<Args>(args)...);

            }
        }
        return try_emplace_base(it, key, std::forward<Args>(args)...);
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

  private:
    /*
     * This method does two things:
     * - It changes the key of the node at 'old_key' to 'new_key'.
     * - It inserts a new node after 'new_key' with value 'old_key'
     * Invariants:
     * - old_key > new_key
     * - there is no other key between new_key and old_key
     */
    NodeT* UpdateKeyAndAddNode(index_t old_key, index_t new_key, NodeT* old_node) {
        assert(old_key > new_key);
        if (data_.size() < M) {
            auto it = lower_bound(old_key);
            assert(it != data_.end());
            it->first = new_key;
            ++it;
            auto new_node = new NodeT(old_node->is_leaf_, this, old_node, old_node->next_node_);
            data_.insert(it, {old_key, new_node});
            if (old_node->next_node_) {
                old_node->next_node_->prev_node_ = new_node;
            }
            return new_node;
        }
        assert(false && "TODO implement");
        return nullptr;
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
            return IterT(this, node, it - node->data_.begin());
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
        return IterT(this, root_);
    }

    [[nodiscard]] auto begin() const {
        return IterT(this, root_);
    }

    [[nodiscard]] auto cbegin() const {
        return IterT(this, root_);
    }

    [[nodiscard]] auto end() {
        return IterT(this, nullptr);
    }

    [[nodiscard]] auto end() const {
        return IterT(this, nullptr);
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
    auto try_emplace_base(index_t key, Args&&... args) {
        auto node = root_;
        while (!node->is_leaf()) {
            auto it = node->lower_bound(key);
            assert(it != node->end());
            node = it->node_;
        }
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

    struct BstStackEntry {
        BstStackEntry(NodeT* node, size_t pos) : node_{node}, pos_{pos} {};

        bool operator==(const BstStackEntry& other) const {
            return node_ == other.node_ && pos_ == other.pos_;
        }

        NodeT* node_;
        size_t pos_;
    };

  public:
    // BstIterator() : first{0}, owner_{nullptr} {};

    // node=nullptr indicates end()
    explicit BstIterator(const b_plus_tree_map<T>* owner, NodeT* node, index_t index = 0)
    : owner_{owner} {
        if (node != nullptr) {
            Push(node, index);
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
            while (!se->node_->is_leaf()) {
                if (se->node_->size() < se->pos_) {
                    ++se->pos_;
                    se = &Push(se->node_->data_[se->pos_].node_);
                }
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

    BstStackEntry& Push(NodeT* node, size_t pos = 0) {
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
        // TODO
        return true;
    }
    const b_plus_tree_map<T>* owner_;
    std::vector<BstStackEntry> stack_{};
    size_t stack_size_{0};
};
}  // namespace

};  // namespace improbable::phtree

#endif  // PHTREE_COMMON_B_PLUS_TREE_H
