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

    enum {
        VALUE = 0,
        NODE = 1,
        EMPTY = 2,
    };

  public:
    /*
     * Construct entry with existing node.
     */
    Entry(const KeyT& k, std::unique_ptr<NodeT>&& node_ptr) noexcept
    : kd_key_{k}
    , node_{std::move(node_ptr)}
    , type{NODE}
#ifdef PH_TREE_ENTRY_POSTLEN
    , postfix_len_{node_->GetPostfixLen()}
#endif
    {
    }

    /*
     * Construct entry with a new node.
     */
    Entry(bit_width_t infix_len, bit_width_t postfix_len) noexcept
    : kd_key_()
    , node_{std::make_unique<NodeT>(infix_len, postfix_len)}
    , type{NODE}
#ifdef PH_TREE_ENTRY_POSTLEN
    , postfix_len_{postfix_len}
#endif
    {
    }

    /*
     * Construct entry with existing T.
     */
    Entry(const KeyT& k, std::optional<ValueT>&& value) noexcept
    : kd_key_{k}
    , value_{std::move(value)}
    , type{VALUE}
#ifdef PH_TREE_ENTRY_POSTLEN
    , postfix_len_{0}
#endif
    {
        //        value.reset();  // std::optional's move constructor does not destruct the previous
    }

    /*
     * Construct entry with new T or moved T.
     */
    template <typename... Args>
    explicit Entry(const KeyT& k, Args&&... args) noexcept
    : kd_key_{k}
    , value_{std::in_place, std::forward<Args>(args)...}
    , type{VALUE}
#ifdef PH_TREE_ENTRY_POSTLEN
    , postfix_len_{0}
#endif
    {
    }

    Entry(const Entry& other) = delete;
    Entry& operator=(const Entry& other) = delete;

    Entry(Entry&& other) noexcept : kd_key_{std::move(other.kd_key_)}, type{std::move(other.type)} {
#ifdef PH_TREE_ENTRY_POSTLEN
        postfix_len_ = std::move(other.postfix_len_);
#endif
        AssignUnion(std::move(other));
    }

    Entry& operator=(Entry&& other) noexcept {
        kd_key_ = std::move(other.kd_key_);
#ifdef PH_TREE_ENTRY_POSTLEN
        postfix_len_ = std::move(other.postfix_len_);
#endif
        DestroyUnion();
        AssignUnion(std::move(other));
        return *this;
    }

    ~Entry() noexcept {
        DestroyUnion();
    }

    [[nodiscard]] const KeyT& GetKey() const {
        return kd_key_;
    }

    [[nodiscard]] bool IsValue() const {
        return type == VALUE;
    }

    [[nodiscard]] bool IsNode() const {
        return type == NODE;
    }

    [[nodiscard]] T& GetValue() const {
        assert(type == VALUE);
        return const_cast<T&>(*value_);
    }

    [[nodiscard]] NodeT& GetNode() const {
        assert(type == NODE);
        return *node_;
    }

    void SetNode(std::unique_ptr<NodeT>&& node) noexcept {
#ifdef PH_TREE_ENTRY_POSTLEN
        postfix_len_ = node->GetPostfixLen();
#endif
        //        std::cout << "size EV : " << sizeof(kd_key_) << " +  " << sizeof(node_) << " +  "
        //                  << sizeof(value_) << "+" << sizeof(type) << " = " << sizeof(*this) <<
        //                  std::endl;
        DestroyUnion();
        type = NODE;
        new (&node_) std::unique_ptr<NodeT>{std::move(node)};
        assert(!node);
    }

    [[nodiscard]] bit_width_t GetNodePostfixLen() const {
        assert(IsNode());
#ifdef PH_TREE_ENTRY_POSTLEN
        return postfix_len_;
#else
        return GetNode().GetPostfixLen();
#endif
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

    void ReplaceNodeWithDataFromEntry(Entry&& other) {
        assert(IsNode());
        // 'other' may be referenced from the local node, so we need to do move(other)
        // before destructing the local node.
        auto node = std::move(node_);
        type = EMPTY;
        *this = std::move(other);
        node.~unique_ptr();
#ifdef PH_TREE_ENTRY_POSTLEN
        postfix_len_ = std::move(other.postfix_len_);
#endif
    }

  private:
    void AssignUnion(Entry&& other) noexcept {
        type = std::move(other.type);
        if (type == NODE) {
            new (&node_) std::unique_ptr<NodeT>{std::move(other.node_)};
        } else if (type == VALUE) {
            new (&value_) std::optional<ValueT>{std::move(other.value_)};
        } else {
            assert(false && "Assigning from an EMPTY variant is a waste of time.");
        }
    }

    void DestroyUnion() noexcept {
        if (type == VALUE) {
            value_.~optional();
        } else if (type == NODE) {
            node_.~unique_ptr();
        } else {
            assert(EMPTY);
        }
        type = EMPTY;
    }

    KeyT kd_key_;
    union {
        std::unique_ptr<NodeT> node_;
        std::optional<ValueT> value_;
    };
    alignas(2) std::uint16_t type;
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

}  // namespace improbable::phtree::v16

#endif  // PHTREE_V16_ENTRY_H
