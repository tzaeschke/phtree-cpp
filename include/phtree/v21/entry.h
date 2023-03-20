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

#ifndef PHTREE_V21_ENTRY_H
#define PHTREE_V21_ENTRY_H

#include "node.h"
#include "phtree/common/common.h"
#include <cassert>
#include <memory>

namespace improbable::phtree::v21 {

template <dimension_t DIM, typename T, typename SCALAR>
class Node;

/*
 * Nodes in the PH-Tree contain up to 2^DIM Entries, one in each geometric quadrant.
 * Entries can contain two types of data:
 * - A key/value pair (value of type T)
 * - A prefix/child-node pair, where prefix is the prefix of the child node and the
 *   child node is contained in a unique_ptr.
 */
template <dimension_t DIM, typename T, typename SCALAR>
class Entry {
    using bit_width_t = detail::bit_width_t;
    using KeyT = PhPoint<DIM, SCALAR>;
    using ValueT = std::remove_const_t<T>;
    using NodeT = Node<DIM, T, SCALAR>;

  public:
    Entry() : node_{nullptr} {}

    /*
     * Construct entry with existing node.
     */
    Entry(const KeyT& k, NodeT* node, bit_width_t postfix_len) noexcept
    : kd_key_{k}
    , node_{node}
    , postfix_len_{static_cast<std::uint16_t>(postfix_len)} {}  // TODO remove 16 bit nonsense

    /*
     * Construct entry with existing T (T is not movable).
     */
//    template <typename ValueT2 = ValueT>
//    Entry(
//        const KeyT& k,
//        ValueT2&& value,
//        typename std::enable_if_t<!std::is_move_constructible_v<ValueT2>, int>::type = 0) noexcept
//    : kd_key_{k}, value_(value), node_{nullptr}, postfix_len_{0} {}

    /*
     * Construct entry with existing T (T must be movable).
     */
    template <typename ValueT2 = ValueT>
    Entry(
        const KeyT& k,
        ValueT2&& value, // TODO backport removed ::type !!!!!
        typename std::enable_if_t<std::is_move_constructible_v<ValueT2>, int> = 0) noexcept
    : kd_key_{k}, value_(std::forward<ValueT2>(value)), node_{nullptr}, postfix_len_{0} {}

    /*
     * Construct entry with new T or copied T (T is not movable).
     */
    template <
        typename ValueT2 = ValueT,
        typename = std::enable_if_t<!std::is_move_constructible_v<ValueT2>>>
    explicit Entry(const KeyT& k, const ValueT& value) noexcept
    : kd_key_{k}, value_(value), node_{nullptr}, postfix_len_{0} {}

    /*
     * Construct entry with new T or copied T (T is not movable, using T's default constructor).
     */
//    template <
//        typename ValueT2 = ValueT,
//        typename = std::enable_if_t<!std::is_move_constructible_v<ValueT2>>>
//    explicit Entry(const KeyT& k) noexcept
//    : kd_key_{k}, value_(), node_{nullptr}, postfix_len_{0} {}

    /*
     * Construct entry with new T or moved T (T must be movable).
     */
//    template <
//        typename... Args,
//        typename ValueT2 = ValueT,
//        typename = std::enable_if_t<std::is_move_constructible_v<ValueT2>>>
//    explicit Entry(const KeyT& k, Args&&... args) noexcept
//    : kd_key_{k}, value_(std::forward<Args>(args)...), node_{nullptr}, postfix_len_{0} {}

    Entry(const Entry& other) = delete;
    Entry& operator=(const Entry& other) = delete;
    Entry(Entry&& other) noexcept :
        kd_key_{std::move(other.kd_key_)},
        node_{other.node_},
        value_{std::move(other.value_)},
        postfix_len_{other.postfix_len_} {
        other.node_ = nullptr;
    }
    Entry& operator=(Entry&& other) noexcept {
        kd_key_ = std::move(other.kd_key_);
        node_ = other.node_;
        value_ = std::move(other.value_);
        postfix_len_ = other.postfix_len_;
        other.node_ = nullptr;
        return *this;
    }
    ~Entry() noexcept {
        if (node_ != nullptr) {
            delete node_;
            node_ = nullptr;
        }
        assert(!is_deleted);
        is_deleted = true;
    }

    bool is_deleted = false;

    void SetNodeCenter() {
        // The node center is defined as the prefix + a '1' bit after the prefix. The remaining
        // bits, i.e. all post_len bits must be '0'.
        // This is required for window queries which would otherwise need to calculate the
        // center each time they traverse a node.
        assert(IsNode());
        detail::bit_mask_t<SCALAR> maskHcBit = detail::bit_mask_t<SCALAR>(1) << postfix_len_;
        detail::bit_mask_t<SCALAR> maskVT = detail::MAX_MASK<SCALAR> << postfix_len_;
        // to prevent problems with signed long when using 64 bit
        if (postfix_len_ < detail::MAX_BIT_WIDTH<SCALAR> - 1) {
            for (dimension_t i = 0; i < DIM; ++i) {
                kd_key_[i] = (kd_key_[i] | maskHcBit) & maskVT;
            }
        } else {
            for (dimension_t i = 0; i < DIM; ++i) {
                kd_key_[i] = 0;
            }
        }
    }

    [[nodiscard]] const KeyT& GetKey() const noexcept {
        return kd_key_;
    }

    [[nodiscard]] bool IsValue() const noexcept {
        return node_ == nullptr;
    }

    [[nodiscard]] bool IsNode() const noexcept {
        return node_ != nullptr;
    }

    [[nodiscard]] bool IsRoot() const noexcept {
        return postfix_len_ == detail::MAX_BIT_WIDTH<SCALAR> - 1;
    }

    [[nodiscard]] T& GetValue() const {
        assert(node_ == nullptr);
        return const_cast<T&>(value_);
    }

    [[nodiscard]] const NodeT& GetNode() const {
        assert(node_ != nullptr);
        return *node_;
    }

    [[nodiscard]] NodeT& GetNode() {
        assert(node_ != nullptr);
        return *node_;
    }

    void SetKey(const KeyT& key) noexcept {
        assert(node_ == nullptr);  // Do we have any other use?
        kd_key_ = key;
    }

    void SetNode(NodeT* node, bit_width_t postfix_len) noexcept {
        postfix_len_ = static_cast<std::uint16_t>(postfix_len);
        // TODO if NODE or VALUE, destruct value?
        assert(node_ == nullptr);
        //        if (node_ != nullptr) {
        //            delete node_;
        //        }
        node_ = node;
        SetNodeCenter();
    }

    [[nodiscard]] bit_width_t GetNodePostfixLen() const noexcept {
        assert(IsNode());
        return postfix_len_;
    }

    [[nodiscard]] bit_width_t GetNodeInfixLen(bit_width_t parent_postfix_len) const noexcept {
        assert(IsNode());
        return parent_postfix_len - GetNodePostfixLen() - 1;
    }

    [[nodiscard]] bool HasNodeInfix(bit_width_t parent_postfix_len) const noexcept {
        assert(IsNode());
        return parent_postfix_len - GetNodePostfixLen() - 1 > 0;
    }

    [[nodiscard]] ValueT&& ExtractValue() noexcept {
        assert(IsValue());
        return std::move(value_);
    }

    [[nodiscard]] NodeT* ExtractNode() noexcept {
        assert(IsNode());
        // Moving the node somewhere else means we should remove it here:
        auto node = node_;
        node_ = nullptr;
        return node;
    }

    void ReplaceNodeWithDataFromEntry(Entry&& other) {
        assert(IsNode());
        // 'other' may be referenced from the local node, so we need to do move(other)
        // before destructing the local node.
        auto node = std::move(node_);
        *this = std::move(other);
        if (IsNode()) {
            SetNodeCenter();
        }
        if (node != nullptr) {
            delete node;
        }
    }

  private:
    KeyT kd_key_;
    NodeT* node_;
    ValueT value_;
    // The length (number of bits) of post fixes (the part of the coordinate that is 'below' the
    // current node). If a variable prefix_len would refer to the number of bits in this node's
    // prefix, and if we assume 64 bit values, the following would always hold:
    // prefix_len + 1 + postfix_len = 64.
    // The '+1' accounts for the 1 bit that is represented by the local node's hypercube,
    // i.e. the same bit that is used to create the lookup keys in entries_.
    std::uint16_t postfix_len_{};
};

}  // namespace improbable::phtree::v21

#endif  // PHTREE_V21_ENTRY_H
