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

#ifndef PHTREE_V16_ITERATOR_FULL_POS_H
#define PHTREE_V16_ITERATOR_FULL_POS_H

#include "iterator_base.h"
#include "node.h"
#include "phtree/common/common.h"

namespace improbable::phtree::v16 {

template <typename T, typename CONVERT>
// TODO FILTER + CONVERTER ????
class IteratorWithFullPos : public IteratorWithFilter<T, CONVERT> {
    static constexpr dimension_t DIM = CONVERT::DimInternal;
    using SCALAR = typename CONVERT::ScalarInternal;
    using EntryT = typename IteratorWithFilter<T, CONVERT>::EntryT;
    using NodeIter = EntryIterator<DIM, EntryT>;
    friend PhTreeV16<DIM, T, CONVERT>;

  public:
    explicit IteratorWithFullPos(const CONVERT* converter) noexcept
    : IteratorWithFilter<T, CONVERT>(nullptr, converter)
    , current_node_{nullptr}
    , parent_node_{nullptr}
    , node_iter_{} {}

    explicit IteratorWithFullPos(
        const EntryT* current_node, const EntryT* parent_node, const CONVERT* converter) noexcept
    : IteratorWithFilter<T, CONVERT>(nullptr, converter)
    , current_node_{current_node}
    , parent_node_{parent_node}
    , node_iter_{} {}

    explicit IteratorWithFullPos(
        NodeIter current_iter,
        const EntryT* current_node,
        const EntryT* parent_node,
        const CONVERT* converter) noexcept
    : IteratorWithFilter<T, CONVERT>(&current_iter->second, converter)
    , current_node_{current_node}
    , parent_node_{parent_node}
    , node_iter_{current_iter} {}

    IteratorWithFullPos& operator++() {
        this->SetFinished();
        return *this;
    }

    IteratorWithFullPos operator++(int) {
        IteratorWithParent iterator(*this);
        ++(*this);
        return iterator;
    }

  private:
    /*
     * The parent entry contains the parent node. The parent node is the node ABOVE the current node
     * which contains the current entry.
     */
    EntryT* GetNodeEntry() const {
        return const_cast<EntryT*>(current_node_);
    }

    EntryT* GetParentNodeEntry() const {
        return const_cast<EntryT*>(parent_node_);
    }

    NodeIter& GetEntryIter() const {
        return const_cast<NodeIter&>(node_iter_);
    }

    const EntryT* current_node_;
    const EntryT* parent_node_;
    NodeIter node_iter_;
};

}  // namespace improbable::phtree::v16

#endif  // PHTREE_V16_ITERATOR_FULL_POS_H
