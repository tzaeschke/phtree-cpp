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

#ifndef PHTREE_COMMON_B_PLUS_TREE_HASH_MAP_H
#define PHTREE_COMMON_B_PLUS_TREE_HASH_MAP_H

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
 * The b_plus_tree_hash_map is a B+tree implementation that uses a hierarchy of horizontally
 * connected nodes for fast traversal through all entries.
 *
 *
 * Behavior:
 * This is a hash map. It uses hashes as key but of course allows multiple entries per key as long
 * as the values are no equal.
 *
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
template <typename T, typename HashT = std::hash<T>, typename PredT = std::equal_to<T>, std::uint64_t COUNT_MAX = 1000000>
class b_plus_tree_hash_set {
    class bpt_node_base;
    template <typename ThisT, typename EntryT>
    class bpt_node_data;
    class bpt_node_leaf;
    class bpt_node_inner;
    class bpt_iterator;

    using hash_t = std::uint64_t;

    using bpt_entry_inner = std::pair<hash_t, bpt_node_base*>;
    using bpt_entry_leaf = std::pair<hash_t, T>;

    using IterT = bpt_iterator;
    using NodeT = bpt_node_base;
    using NLeafT = bpt_node_leaf;
    using NInnerT = bpt_node_inner;
    using LeafIteratorT = decltype(std::vector<bpt_entry_leaf>().begin());
    using TreeT = b_plus_tree_hash_set<T, HashT, PredT, COUNT_MAX>;

  public:
    explicit b_plus_tree_hash_set() : root_{new NLeafT(nullptr, nullptr, nullptr)}, size_{0} {};

    b_plus_tree_hash_set(const b_plus_tree_hash_set& other) : size_{other.size_} {
        root_ = other.root_->is_leaf() ? new NLeafT(*other.root_->as_leaf())
                                       : new NInnerT(*other.root_->as_inner());
    }

    b_plus_tree_hash_set(b_plus_tree_hash_set&& other) noexcept
    : root_{other.root_}, size_{other.size_} {
        other.root_ = nullptr;
        other.size_ = 0;
    }

    b_plus_tree_hash_set& operator=(const b_plus_tree_hash_set& other) {
        assert(this != &other);
        delete root_;
        root_ = other.root_->is_leaf() ? new NLeafT(*other.root_->as_leaf())
                                       : new NInnerT(*other.root_->as_inner());
        size_ = other.size_;
        return *this;
    }

    b_plus_tree_hash_set& operator=(b_plus_tree_hash_set&& other) noexcept {
        delete root_;
        root_ = other.root_;
        other.root_ = nullptr;
        size_ = other.size_;
        other.size_ = 0;
        return *this;
    }

    ~b_plus_tree_hash_set() {
        delete root_;
        root_ = nullptr;
    }

    [[nodiscard]] auto find(const T& value) noexcept {
        auto node = root_;
        auto hash = HashT{}(value);
        while (!node->is_leaf()) {
            node = node->as_inner()->find(hash);
            if (node == nullptr) {
                return end();
            }
        }
        return node->as_leaf()->find(hash, value);
    }

    [[nodiscard]] auto find(const T& value) const noexcept {
        return const_cast<b_plus_tree_hash_set&>(*this).find(value);
    }

    [[nodiscard]] auto lower_bound(const T& value) noexcept {
        auto node = root_;
        auto hash = HashT{}(value);
        while (!node->is_leaf()) {
            node = node->as_inner()->find(hash);
            if (node == nullptr) {
                return end();
            }
        }
        return node->as_leaf()->lower_bound_as_iter(hash, value);
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
        T t(std::forward<Args>(args)...);
        hash_t hash = HashT{}(t);
        return try_emplace_base(hash, std::move(t));
    }

    template <typename... Args>
    auto emplace_hint(const IterT& hint, Args&&... args) {
        T t(std::forward<Args>(args)...);
        auto hash = HashT{}(t);
        return try_emplace_base_hint(hint, hash, std::move(t));
    }

    size_t erase(const T& value) {
        auto node = root_;
        auto hash = HashT{}(value);
        while (!node->is_leaf()) {
            node = node->as_inner()->find(hash);
            if (node == nullptr) {
                return 0;
            }
        }
        auto n = node->as_leaf()->erase_key(hash, value, *this);
        size_ -= n;
        return n;
    }

    void erase(const IterT& iterator) {
        assert(iterator != end());
        --size_;
        iterator.node_->erase_it(iterator.iter_, *this);
    }

    [[nodiscard]] size_t size() const noexcept {
        return size_;
    }

    void _check() {
        size_t count = 0;
        NLeafT* prev_leaf = nullptr;
        hash_t known_min = std::numeric_limits<hash_t>::max();
        root_->_check(count, nullptr, prev_leaf, known_min, 0);
        assert(count == size());
    }

  private:
    template <typename T2 = T>
    auto try_emplace_base(hash_t hash, T2&& t) {
        auto node = root_;
        while (!node->is_leaf()) {
            node = node->as_inner()->find_or_last(hash);
        }
        return node->as_leaf()->try_emplace(hash, *this, size_, std::forward<T2>(t));
    }

    template <typename T2> // TODO universal ref?? -> Probably no, we always pass in a rvalue...
    auto try_emplace_base_hint(const IterT& hint, hash_t hash, T2&& t) {
        // TODO sanity checkes: correct node?    ==end()?
        // TODO use hint inside node
        // TODO can we use this in the PH-TreeV16?
        auto node = hint.node_;
        return node->as_leaf()->try_emplace(hash, *this, size_, std::forward<T2>(t));
    }

    class bpt_node_base {
      public:
        explicit bpt_node_base(bool is_leaf, NInnerT* parent) noexcept
        : is_leaf_{is_leaf}, parent_{parent} {}

        virtual ~bpt_node_base() noexcept = default;

        [[nodiscard]] inline bool is_leaf() const noexcept {
            return is_leaf_;
        }

        [[nodiscard]] inline NInnerT* as_inner() noexcept {
            assert(!is_leaf_);
            return static_cast<NInnerT*>(this);
        }

        [[nodiscard]] inline NLeafT* as_leaf() noexcept {
            assert(is_leaf_);
            return static_cast<NLeafT*>(this);
        }

        virtual void _check(size_t&, NInnerT*, NLeafT*&, hash_t&, hash_t) = 0;

      public:
        const bool is_leaf_;
        NInnerT* parent_;
    };

    template <typename ThisT, typename EntryT>
    class bpt_node_data : public bpt_node_base {
        using DataIteratorT = decltype(std::vector<EntryT>().begin());
        friend IterT;

        constexpr static size_t M_leaf = std::min(size_t(16), COUNT_MAX);
        // Default MAX is 32. Special case for small COUNT with smaller inner leaf or
        // trees with a single inner leaf. '*2' is added because leaf filling is not compact.
        constexpr static size_t M_inner = std::min(size_t(16), COUNT_MAX / M_leaf * 2);
        // TODO This could be improved but requires a code change to move > 1 entry when merging.
        constexpr static size_t M_leaf_min = 2;   // std::max((size_t)2, M_leaf >> 2);
        constexpr static size_t M_inner_min = 2;  // std::max((size_t)2, M_inner >> 2);
        // There is no point in allocating more leaf space than the max amount of entries.
        constexpr static size_t M_leaf_init = std::min(size_t(8), COUNT_MAX);
        constexpr static size_t M_inner_init = 4;

      public:
        explicit bpt_node_data(bool is_leaf, NInnerT* parent, ThisT* prev, ThisT* next) noexcept
        : bpt_node_base(is_leaf, parent), data_{}, prev_node_{prev}, next_node_{next} {
            data_.reserve(this->M_init());
        }

        virtual ~bpt_node_data() noexcept = default;

        [[nodiscard]] inline size_t M_min() {
            return this->is_leaf_ ? M_leaf_min : M_inner_min;
        }

        [[nodiscard]] inline size_t M_max() {
            return this->is_leaf_ ? M_leaf : M_inner;
        }

        [[nodiscard]] inline size_t M_init() {
            return this->is_leaf_ ? M_leaf_init : M_inner_init;
        }

        [[nodiscard]] auto lower_bound(hash_t hash) noexcept {
            return std::lower_bound(
                data_.begin(), data_.end(), hash, [](EntryT& left, const hash_t hash) {
                    return left.first < hash;
                });
        }

        [[nodiscard]] size_t size() const noexcept {
            return data_.size();
        }

        void erase_entry(DataIteratorT it_to_erase, TreeT& tree) {
            auto& parent_ = this->parent_;
            hash_t max_key_old = data_.back().first;

            data_.erase(it_to_erase);
            if (parent_ == nullptr) {
                if constexpr (std::is_same_v<ThisT, NInnerT>) {
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

            if (data_.empty()) {
                // Nothing to merge, just remove node. This should be rare, i.e. only happens when
                // a rare 1-entry node has its last entry removed.
                remove_from_siblings();
                parent_->remove_node(max_key_old, this, tree);
                return;
            }

            if (data_.size() < this->M_min()) {
                // merge
                if (prev_node_ != nullptr && prev_node_->data_.size() < this->M_max()) {
                    remove_from_siblings();
                    auto& prev_data = prev_node_->data_;
                    if constexpr (std::is_same_v<ThisT, NLeafT>) {
                        prev_data.emplace_back(std::move(data_[0]));
                    } else {
                        data_[0].second->parent_ = prev_node_;
                        prev_data.emplace_back(std::move(data_[0]));
                        data_[0].second = nullptr;
                    }
                    auto prev_node = prev_node_;  // create copy because (this) will be deleted
                    parent_->remove_node(max_key_old, this, tree);
                    if (prev_node->parent_ != nullptr) {
                        hash_t old1 = (prev_data.end() - 2)->first;
                        hash_t new1 = (prev_data.end() - 1)->first;
                        prev_node->parent_->update_key(old1, new1);
                    }
                    return;
                } else if (next_node_ != nullptr && next_node_->data_.size() < this->M_max()) {
                    remove_from_siblings();
                    auto& next_data = next_node_->data_;
                    if constexpr (std::is_same_v<ThisT, NLeafT>) {
                        next_data.emplace(next_data.begin(), std::move(data_[0]));
                    } else {
                        data_[0].second->parent_ = next_node_;
                        next_data.emplace(next_data.begin(), std::move(data_[0]));
                        data_[0].second = nullptr;
                    }
                    parent_->remove_node(max_key_old, this, tree);
                    return;
                }
                // This node is too small but there is nothing we can do.
            }
            if (it_to_erase == data_.end()) {
                parent_->update_key(max_key_old, data_.back().first);
            }
        }

        auto prepare_emplace(hash_t key, TreeT& tree, DataIteratorT& it_in_out) {
            if (data_.size() < this->M_max()) {
                if (this->parent_ != nullptr && key > data_.back().first) {
                    this->parent_->update_key(data_.back().first, key);
                }
                return static_cast<ThisT*>(this);
            }

            ThisT* dest = this->split_node(key, tree);
            if (dest != this) {
                // The insertion pos in node2 can be calculated:
                auto old_pos = it_in_out - data_.begin();
                it_in_out = dest->data_.begin() + old_pos - data_.size();
            }
            return dest;
        }

        void _check_data(NInnerT* parent, hash_t known_max) {
            (void)parent;
            (void)known_max;
            // assert(parent_ == nullptr || data_.size() >= M_min);
            assert(this->parent_ == parent);
            if (this->data_.empty()) {
                assert(parent == nullptr);
                return;
            }
            assert(this->parent_ == nullptr || known_max == this->data_.back().first);
        }

      private:
        ThisT* split_node(hash_t key, TreeT& tree) {
            auto max_key = data_.back().first;
            if (this->parent_ == nullptr) {
                auto* new_parent = new NInnerT(nullptr, nullptr, nullptr);
                new_parent->emplace_back(max_key, this);
                tree.root_ = new_parent;
                this->parent_ = new_parent;
            }

            // create new node
            auto* node2 = new ThisT(this->parent_, static_cast<ThisT*>(this), next_node_);
            if (next_node_ != nullptr) {
                next_node_->prev_node_ = node2;
            }
            next_node_ = node2;

            // populate new node
            // TODO Optimize populating new node: move 1st part, insert new value, move 2nd part...?
            auto split_pos = this->M_max() >> 1;
            node2->data_.insert(
                node2->data_.end(),
                std::make_move_iterator(data_.begin() + split_pos),
                std::make_move_iterator(data_.end()));
            data_.erase(data_.begin() + split_pos, data_.end());

            if constexpr (std::is_same_v<ThisT, NInnerT>) {
                for (auto& e : node2->data_) {
                    e.second->parent_ = node2;
                }
            }

            // Add node to parent
            auto split_key = data_[split_pos - 1].first;
            if (key > split_key && key < node2->data_[0].first) {
                // This is a bit hacky:
                // Add new entry at END of first node when possible -> avoids some shifting
                split_key = key;
            }
            this->parent_->update_key_and_add_node(
                max_key, split_key, std::max(max_key, key), node2, tree);

            // Return node for insertion of new value
            return key > split_key ? node2 : static_cast<ThisT*>(this);
        }

        void remove_from_siblings() {
            if (next_node_ != nullptr) {
                next_node_->prev_node_ = prev_node_;
            }
            if (prev_node_ != nullptr) {
                prev_node_->next_node_ = next_node_;
            }
        }

      protected:
        std::vector<EntryT> data_;
        ThisT* prev_node_;
        ThisT* next_node_;
    };

    class bpt_node_leaf : public bpt_node_data<bpt_node_leaf, bpt_entry_leaf> {
      public:
        explicit bpt_node_leaf(NInnerT* parent, NLeafT* prev, NLeafT* next) noexcept
        : bpt_node_data<NLeafT, bpt_entry_leaf>(true, parent, prev, next) {}

        ~bpt_node_leaf() noexcept = default;

//        [[nodiscard]] IterT find(hash_t hash) noexcept {
//            auto it = this->lower_bound(hash);
//            if (it != this->data_.end() && it->first == key) {
//                return IterT(this, it);
//            }
//            return IterT();
//        }

        [[nodiscard]] IterT find(hash_t hash, const T& value) noexcept {
            auto it = this->lower_bound(hash);
            PredT equals{};
            IterT iter_full(this, it);
            while (!iter_full.is_end() && iter_full.hash() == hash) {
                if (equals(*iter_full, value)) {
                    return iter_full;
                }
                ++iter_full;
            }
            return IterT();
        }

        [[nodiscard]] IterT lower_bound_as_iter(hash_t hash) noexcept {
            auto it = this->lower_bound(hash);
            if (it != this->data_.end()) {
                return IterT(this, it);
            }
            return IterT();
        }

        // TODO remobve/rename !?!?!?
        [[nodiscard]] auto lower_bound2(hash_t hash, const T& value) noexcept {
            auto it = this->lower_bound(hash);
            IterT iter_full(this, it);
            while (!iter_full.is_end() && iter_full.hash() == hash) {
                if (PredT{}(*iter_full, value)) {
                    break;
                }
                ++iter_full;
            }
            return iter_full;
        }

        template <typename... Args>
        auto try_emplace(hash_t hash, TreeT& tree, size_t& entry_count, Args&&... args) {
            auto it = this->lower_bound(hash);

            if (it == this->data_.end() || it->first != hash) {
                // fast track insert
                ++entry_count;
                auto dest = this->prepare_emplace(hash, tree, it);
                auto x = dest->data_.emplace(
                    it,
                    std::piecewise_construct,
                    std::forward_as_tuple(hash),
                    std::forward_as_tuple(std::forward<Args>(args)...));
                return std::make_pair(IterT(this, x), true);
            }

            // Hash collision !
            T t(std::forward<Args>(args)...);
            PredT equals{}; // static?

            IterT full_iter(this, it);
            while (!full_iter.is_end() && full_iter.hash() == hash) {
                if (equals(*full_iter, t)) {
                    return std::make_pair(full_iter, false);
                }
                ++full_iter;
            }
//            if (it != this->data_.end() && it->first == key) {
//                return std::make_pair(IterT(this, it), false);
//            }
            ++entry_count;

            auto dest = this->prepare_emplace(hash, tree, it);

            //            auto x = dest->data_.emplace(
            //                it,
            //                std::piecewise_construct,
            //                std::forward_as_tuple(key),
            //                std::forward_as_tuple(std::forward<Args>(args)...));
            auto x = dest->data_.emplace(
                it,  // TODO "it"? Should we use sorted insert?
                std::piecewise_construct,
                std::forward_as_tuple(hash),
                std::forward_as_tuple(std::move(t)));
            return std::make_pair(IterT(this, x), true);
        }

        bool erase_key(hash_t hash, const T& value, TreeT& tree) {
            auto iter_node = this->lower_bound(hash);
            IterT iter_full(this, iter_node);
            while (!iter_full.is_end() && iter_full->first == hash) {
                if (PredT{}(*iter_full, value)) {
                    this->erase_entry(iter_full.iter_, tree);
                    return true;
                }
            }
            return false;
        }

        void erase_it(LeafIteratorT iter, TreeT& tree) {
            this->erase_entry(iter, tree);
        }

        void _check(
            size_t& count, NInnerT* parent, NLeafT*& prev_leaf, hash_t& known_min, hash_t known_max) {
            this->_check_data(parent, known_max);

            assert(prev_leaf == this->prev_node_);
            for (auto& e : this->data_) {
                assert(count == 0 || e.first >= known_min);
                assert(this->parent_ == nullptr || e.first <= known_max);
                ++count;
                known_min = e.first;
            }
            prev_leaf = this;
        }
    };

    class bpt_node_inner : public bpt_node_data<bpt_node_inner, bpt_entry_inner> {
      public:
        explicit bpt_node_inner(NInnerT* parent, NInnerT* prev, NInnerT* next) noexcept
        : bpt_node_data<NInnerT, bpt_entry_inner>(false, parent, prev, next) {}

        ~bpt_node_inner() noexcept {
            for (auto& e : this->data_) {
                if (e.second != nullptr) {
                    delete e.second;
                }
            }
        }

        [[nodiscard]] NodeT* find(hash_t hash) noexcept {
            auto it = this->lower_bound(hash);
            return it != this->data_.end() ? it->second : nullptr;
        }

        [[nodiscard]] NodeT* find_or_last(hash_t hash) noexcept {
            auto it = this->lower_bound(hash);
            return it != this->data_.end() ? it->second : this->data_.back().second;
        }

        void emplace_back(hash_t hash, NodeT* node) {
            this->data_.emplace_back(hash, node);
        }

        void _check(
            size_t& count, NInnerT* parent, NLeafT*& prev_leaf, hash_t& known_min, hash_t known_max) {
            this->_check_data(parent, known_max);

            assert(this->parent_ == nullptr || known_max == this->data_.back().first);
            auto prev_key = this->data_[0].first;
            int n = 0;
            for (auto& e : this->data_) {
                assert(n == 0 || e.first > prev_key);
                e.second->_check(count, this, prev_leaf, known_min, e.first);
                assert(this->parent_ == nullptr || e.first <= known_max);
                prev_key = e.first;
                ++n;
            }
        }

        void update_key(hash_t old_key, hash_t new_key) {
            assert(new_key != old_key);
            auto it = this->lower_bound(old_key);
            assert(it != this->data_.end());
            assert(it->first == old_key);
            it->first = new_key;
            if (this->parent_ != nullptr && ++it == this->data_.end()) {
                this->parent_->update_key(old_key, new_key);
            }
        }

        /*
         * This method does two things:
         * - It changes the key of the node (node 1) at 'key1_old' to 'key1_new'.
         * - It inserts a new node (node 2) after 'new_key1' with value 'key2'
         * Invariants:
         * - Node1: key1_old > key1_new; Node 1 vs 2: key2 > new_key1
         */
        void update_key_and_add_node(
            hash_t key1_old, hash_t key1_new, hash_t key2, NodeT* child2, TreeT& tree) {
            assert(key2 > key1_new);
            assert(key1_old >= key1_new);
            auto it2 = this->lower_bound(key1_old) + 1;

            auto dest = this->prepare_emplace(key2, tree, it2);
            // prepare_emplace() guarantees that child2 is in the same node as child1
            assert(it2 != dest->data_.begin());
            (it2 - 1)->first = key1_new;
            child2->parent_ = dest;
            dest->data_.emplace(it2, key2, child2);
        }

        void remove_node(hash_t key_remove, NodeT* node, TreeT& tree) {
            auto it_to_erase = this->lower_bound(key_remove);
            while (it_to_erase != this->data_.end() && it_to_erase->first == key_remove) {
                if (it_to_erase->second == node) {
                    delete it_to_erase->second;
                    this->erase_entry(it_to_erase, tree);
                }
                ++it_to_erase;
            }
            assert(false && "Node not found!");
        }
    };

    class bpt_iterator {
        using EntryT = typename b_plus_tree_hash_set<T, HashT, PredT, COUNT_MAX>::bpt_entry_leaf;
        friend b_plus_tree_hash_set<T, HashT, PredT, COUNT_MAX>;

      public:
        using iterator_category = std::forward_iterator_tag;
        using value_type = T;
        using difference_type = std::ptrdiff_t;
        using pointer = T*;
        using reference = T&;

        // Arbitrary position iterator
        explicit bpt_iterator(NLeafT* node, LeafIteratorT it) noexcept : node_{node}, iter_{it} {
            assert(node->is_leaf_ && "just for consistency, insist that we iterate leaves only ");
        }

        // begin() iterator
        explicit bpt_iterator(NodeT* node) noexcept {
            assert(node->parent_ == nullptr && "must start with root node");
            // move iterator to first value
            while (!node->is_leaf_) {
                node = node->as_inner()->data_[0].second;
            }
            node_ = node->as_leaf();

            if (node_->size() == 0) {
                node_ = nullptr;
                iter_ = {};
                return;
            }
            iter_ = node_->data_.begin();
        }

        // end() iterator
        bpt_iterator() noexcept : node_{nullptr}, iter_{} {}

        auto& operator*() const noexcept {
            assert(AssertNotEnd());
            return const_cast<T&>(iter_->second);
        }

        auto* operator->() const noexcept {
            assert(AssertNotEnd());
            return const_cast<T*>(&iter_->second);
        }

        auto& operator++() noexcept {
            assert(AssertNotEnd());
            ++iter_;
            if (iter_ == node_->data_.end()) {
                // this may be a nullptr -> end of data
                node_ = node_->next_node_;
                iter_ = node_ != nullptr ? node_->data_.begin() : LeafIteratorT{};
            }
            return *this;
        }

        auto operator++(int) const noexcept {
            IterT iterator(*this);
            ++(*this);
            return iterator;
        }

        friend bool operator==(const IterT& left, const IterT& right) noexcept {
            return left.iter_ == right.iter_ && left.node_ == right.node_;
        }

        friend bool operator!=(const IterT& left, const IterT& right) noexcept {
            return !(left == right);
        }

        // TODO private
        bool is_end() const noexcept {
            return node_ == nullptr;
        }
      private:
        [[nodiscard]] inline bool AssertNotEnd() const noexcept {
            return node_ != nullptr;
        }

        hash_t hash() {
            return iter_->first;
        }

        NLeafT* node_;
        LeafIteratorT iter_;
    };

  private:
    NodeT* root_;
    size_t size_;
};

template <typename KeyT, typename ValueT, typename HashT = std::hash<KeyT>, typename PredT = std::equal_to<KeyT>, size_t COUNT_MAX = 10000>
class b_plus_tree_hash_map {
    class iterator;
    using IterT = iterator;
    using EntryT = std::pair<KeyT, ValueT>;

  public:
    b_plus_tree_hash_map() : map_{} {};

    b_plus_tree_hash_map(const b_plus_tree_hash_map&) = default;
    b_plus_tree_hash_map(b_plus_tree_hash_map&&) noexcept = default;
    b_plus_tree_hash_map& operator=(const b_plus_tree_hash_map&) = default;
    b_plus_tree_hash_map& operator=(b_plus_tree_hash_map&&) noexcept = default;
    ~b_plus_tree_hash_map() = default;

    auto begin() const {
        return IterT(map_.begin());
    }

    auto end() const {
        return IterT(map_.end());
    }

    auto find(const KeyT& key) const {
        return IterT(map_.find(EntryT{key, {}}));
    }

    template <typename... Args>
    auto emplace(Args&&... args) {
        return try_emplace(std::forward<Args>(args)...);
    }

    template <typename... Args>
    auto emplace_hint(const IterT& hint, Args&&... args) {
        EntryT key(std::forward<Args>(args)...);
        auto result = map_.emplace_hint(hint.map_iter_, std::move(key)); // TODO forward?
        return iterator(result.first);
    }

    template <typename... Args>
    auto try_emplace(Args&&... args) {
        //EntryT key(std::forward<Args>(args)...);
        //auto result = map_.emplace(std::move(key));// TODO forward?
        auto result = map_.emplace(std::make_pair(std::forward<Args>(args)...));// TODO forward?
        return std::make_pair(iterator(result.first), result.second);
    }

    template <typename... Args>
    auto try_emplace(const IterT& hint, Args&&... args) {
        //EntryT key(std::forward<Args>(args)...);
        EntryT key = std::make_pair(std::forward<Args>(args)...);
        auto result = map_.emplace_hint(hint.map_iter_, std::move(key)); // TODO forward?
        return iterator(result.first);
    }

    auto erase(const KeyT& key) {
        return map_.erase({key, {}});
    }

    auto erase(const IterT& iterator) {
        map_.erase(iterator.map_iter_);
    }

    auto size() const {
        return map_.size();
    }

    auto empty() const {
        return map_.size() == 0;
    }

    void _check() {
        // TODO
    }

  private:
    struct EntryHashT {
        size_t operator()(const EntryT& x) const {
            return HashT{}(x.first);
        }
    };

    struct EntryEqualsT {
        bool operator()(const EntryT& x, const EntryT& y) const {
            return PredT{}(x.first, y.first);
        }
    };

    class iterator {
        using T = EntryT;
        using MapIterType = decltype(std::declval<b_plus_tree_hash_set<EntryT, EntryHashT, EntryEqualsT, COUNT_MAX>>().begin());
        friend b_plus_tree_hash_map<KeyT, ValueT, HashT, PredT, COUNT_MAX>;

      public:
        using iterator_category = std::forward_iterator_tag;
        using value_type = T;
        using difference_type = std::ptrdiff_t;
        using pointer = T*;
        using reference = T&;

        explicit iterator(MapIterType map_iter) noexcept : map_iter_{map_iter} {}

        // end() iterator
        iterator() noexcept : map_iter_{} {}

        auto& operator*() const noexcept {
            return const_cast<T&>(*map_iter_);
        }

        auto* operator->() const noexcept {
            return const_cast<T*>(&*map_iter_);  /// TODO remove cast?
        }

        auto& operator++() noexcept {
            ++map_iter_;
            return *this;
        }

        auto operator++(int) noexcept {
            IterT iterator(*this);
            ++(*this);
            return iterator;
        }

        friend bool operator==(const IterT& left, const IterT& right) noexcept {
            return left.map_iter_ == right.map_iter_;
        }

        friend bool operator!=(const IterT& left, const IterT& right) noexcept {
            return !(left == right);
        }

      private:
        MapIterType map_iter_;
    };

    b_plus_tree_hash_set<EntryT, EntryHashT, EntryEqualsT, COUNT_MAX> map_;
};

}  // namespace improbable::phtree

#endif  // PHTREE_COMMON_B_PLUS_TREE_HASH_MAP_H
