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

#ifndef PHTREE_V16_ENTRY_H
#define PHTREE_V16_ENTRY_H

#include "../../phtree/common/common.h"
#include "node.h"
#include <cassert>
#include <memory>
#include <optional>

//#define PH_TREE_ENTRY_POSTLEN 1

namespace improbable::phtree::v16 {

template <dimension_t DIM, typename T, typename SCALAR>
class Node;

template <dimension_t DIM, typename T, typename SCALAR>
struct EntryVariant;

/*
 * Nodes in the PH-Tree contain up to 2^DIM Entries, one in each geometric quadrant.
 * Entries can contain two types of data:
 * - A key/value pair (value of type T)
 * - A prefix/child-node pair, where prefix is the prefix of the child node and the
 *   child node is contained in a unique_ptr.
 */
template <dimension_t DIM, typename T, typename SCALAR>
class Entry {
    using KeyT = PhPoint<DIM, SCALAR>;
    using ValueT = std::remove_const_t<T>;
    using NodeT = Node<DIM, T, SCALAR>;

  public:
    /*
     * Construct entry with existing node.
     */
    Entry(const KeyT& k, std::unique_ptr<NodeT>&& node_ptr)
    : kd_key_{k}
    , variant_{std::move(node_ptr)}
#ifdef PH_TREE_ENTRY_POSTLEN
    , postfix_len_{variant_.GetNode().GetPostfixLen()}
#endif
    {
    }

    /*
     * Construct entry with a new node.
     */
    Entry(bit_width_t infix_len, bit_width_t postfix_len)
    : kd_key_()
    , variant_{std::make_unique<NodeT>(infix_len, postfix_len)}
#ifdef PH_TREE_ENTRY_POSTLEN
    , postfix_len_{postfix_len}
#endif
    {
    }

    /*
     * Construct entry with existing T.
     */
    Entry(const KeyT& k, std::optional<ValueT>&& value)
    : kd_key_{k}
    , variant_{std::move(value)}
#ifdef PH_TREE_ENTRY_POSTLEN
    , postfix_len_{0}
#endif
    {
    }

    /*
     * Construct entry with new T or moved T.
     */
    template <typename... Args>
    explicit Entry(const KeyT& k, Args&&... args)
    : kd_key_{k}
    , variant_{std::in_place, std::forward<Args>(args)...}
#ifdef PH_TREE_ENTRY_POSTLEN
    , postfix_len_{0}
#endif
    {
    }

    [[nodiscard]] const KeyT& GetKey() const {
        return kd_key_;
    }

    [[nodiscard]] bool IsValue() const {
        return variant_.IsValue();
    }

    [[nodiscard]] bool IsNode() const {
        return variant_.IsNode();
    }

    [[nodiscard]] T& GetValue() const {
        assert(IsValue());
        return const_cast<T&>(variant_.GetValue());
    }

    [[nodiscard]] NodeT& GetNode() const {
        assert(IsNode());
        return variant_.GetNode();
    }

    void SetNode(std::unique_ptr<NodeT>&& node) {
#ifdef PH_TREE_ENTRY_POSTLEN
        postfix_len_ = node->GetPostfixLen();
#endif
        variant_.SetNode(std::move(node));
    }

    [[nodiscard]] bit_width_t GetNodePostfixLen() const {
        assert(IsNode());
#ifdef PH_TREE_ENTRY_POSTLEN
        return postfix_len_;
#else
        return variant_.GetNode().GetPostfixLen();
#endif
    }

    [[nodiscard]] std::optional<ValueT>&& ExtractValue() {
        assert(IsValue());
        return variant_.ExtractValue();
    }

    [[nodiscard]] std::unique_ptr<NodeT>&& ExtractNode() {
        assert(IsNode());
        return variant_.ExtractNode();
    }

    void ReplaceNodeWithDataFromEntry(Entry&& other) {
        assert(IsNode());
        kd_key_ = other.GetKey();
        variant_ = std::move(other.variant_);
#ifdef PH_TREE_ENTRY_POSTLEN
        if (IsNode()) {
            postfix_len_ = other.postfix_len_;
        }
#endif
    }

  private:
    KeyT kd_key_;
    EntryVariant<DIM, T, SCALAR> variant_;
    // The length (number of bits) of post fixes (the part of the coordinate that is 'below' the
    // current node). If a variable prefix_len would refer to the number of bits in this node's
    // prefix, and if we assume 64 bit values, the following would always hold:
    // prefix_len + 1 + postfix_len = 64.
    // The '+1' accounts for the 1 bit that is represented by the local node's hypercube,
    // i.e. the same bit that is used to create the lookup keys in entries_.
#ifdef PH_TREE_ENTRY_POSTLEN
    alignas(2) bit_width_t postfix_len_;
#endif
};

//#pragma pack(push, 4)
template <dimension_t DIM, typename T, typename SCALAR>
struct EntryVariant {
    using ValueT = std::remove_const_t<T>;
    using NodeT = Node<DIM, T, SCALAR>;

    enum {
        VALUE = 0,
        NODE = 1,
        EMPTY = 2,
    };

    EntryVariant() = delete;

    explicit EntryVariant(std::unique_ptr<NodeT>&& node_ptr)
    : node_{std::move(node_ptr)}, type{NODE} {}

    explicit EntryVariant(std::optional<ValueT>&& value) : value_{std::move(value)}, type{VALUE} {
        assert(IsValue());
        value.reset();  // TODO why???
        assert(!value.has_value());
    }

    template <typename... Args>
    explicit EntryVariant(std::in_place_t, Args&&... args)
    : value_{std::in_place, std::forward<Args>(args)...}, type{VALUE} {
        assert(IsValue());
    }

    EntryVariant(const EntryVariant& other) = delete;
    EntryVariant& operator=(const EntryVariant& other) = delete;

    EntryVariant(EntryVariant&& other) noexcept : type{std::move(other.type)} {
        switch (type) {
        case VALUE:
            new (&value_) std::optional<ValueT>{std::move(other.value_)};
            other.value_.reset();  // TODO why???
            assert(!other.value_.has_value());
            break;
        case NODE:
            new (&node_) std::unique_ptr<NodeT>{std::move(other.node_)};
            break;
        default:
            assert(false);
        }
        other.type = EMPTY;
        assert(type != EMPTY);
    }

    EntryVariant& operator=(EntryVariant&& other) noexcept {
        switch (other.type) {
        case VALUE:
            if (IsNode()) {
                auto node = std::move(node_);
                new (&value_) std::optional<ValueT>{std::move(other.value_)};
                other.value_.reset();  // TODO why???
                // other.value_.~optional();
                other.type = EMPTY;
                assert(!other.value_.has_value());
                node.~unique_ptr();
            } else if (IsValue()) {
                value_ = std::move(other.value_);
            } else {
                new (&value_) std::optional<ValueT>{std::move(other.value_)};
                other.value_.reset();  // TODO why???
                // other.value_.~optional();
                assert(!other.value_.has_value());
                other.type = EMPTY;
            }
            type = VALUE;
            break;
        case NODE:
            if (IsNode()) {
                auto node = std::move(node_);
                node_ = std::move(other.node_);
                //new (&node_) std::unique_ptr<NodeT>{std::move(other.node_)};
                other.type = EMPTY;
                node.~unique_ptr();
            } else {
                Destroy();
                new (&node_) std::unique_ptr<NodeT>{std::move(other.node_)};
                other.type = EMPTY;
            }
            type = NODE;
            break;
        default:
            assert(false);
        }
        assert(type != EMPTY);
        // other.type = EMPTY;
        return *this;
    }

    ~EntryVariant() noexcept {
        Destroy();
    }

    void Destroy() noexcept {
        switch (type) {
        case VALUE:
            value_.~optional();
            break;
        case NODE:
            node_.~unique_ptr();
            break;
        case EMPTY:
            break;
        default:
            assert(false);
        }
        type = EMPTY;
    }

    void SetValue(std::optional<ValueT>&& value) noexcept {
        Destroy();
        type = VALUE;
        new (&value_) std::optional<ValueT>{std::move(value)};
        assert(!value);
    }

    void SetNode(std::unique_ptr<NodeT>&& node) noexcept {
        Destroy();
        type = NODE;
        // move unique_ptr without calling destructor on previous value
        new (&node_) std::unique_ptr<NodeT>{std::move(node)};
        assert(!node);
    }

    [[nodiscard]] bool IsNode() const {
        return type == NODE;
    }

    [[nodiscard]] bool IsValue() const {
        return type == VALUE;
    }

    [[nodiscard]] NodeT& GetNode() const {
        assert(type == NODE);
        return *node_;
    }

    [[nodiscard]] T& GetValue() const {
        assert(type == VALUE);
        return const_cast<T&>(*value_);
    }

    [[nodiscard]] std::optional<ValueT>&& ExtractValue() {
        assert(IsValue());
        type = EMPTY;
        return std::move(value_);
    }

    [[nodiscard]] std::unique_ptr<NodeT>&& ExtractNode() {
        assert(IsNode());
        type = EMPTY;
        return std::move(node_);
    }

  private:
    union {
        std::unique_ptr<NodeT> node_;
        std::optional<ValueT> value_;
    };
    alignas(2) std::uint16_t type;
};
//#pragma pack(pop)

}  // namespace improbable::phtree::v16

#endif  // PHTREE_V16_ENTRY_H
