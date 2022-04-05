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

#ifndef PHTREE_V16_ITERATOR_SIMPLE_H
#define PHTREE_V16_ITERATOR_SIMPLE_H

#include "../common/common.h"
#include "iterator_base.h"

namespace improbable::phtree::v16 {

template <typename T, typename CONVERT>
class IteratorSimple : public IteratorBase<T, CONVERT> {
    static constexpr dimension_t DIM = CONVERT::DimInternal;
    using SCALAR = typename CONVERT::ScalarInternal;
    using EntryT = typename IteratorBase<T, CONVERT>::EntryT;

  public:
    explicit IteratorSimple(const CONVERT* converter) : IteratorBase<T, CONVERT>(converter) {
        this->SetFinished();
    }

    explicit IteratorSimple(
        const EntryT* current_result,
        const EntryT* current_node,
        const EntryT* parent_node,
        const CONVERT* converter)
    : IteratorBase<T, CONVERT>(converter) {
        if (current_result) {
            this->SetCurrentResult(current_result);
            this->SetCurrentNodeEntry(current_node);
            this->SetParentNodeEntry(parent_node);
        } else {
            this->SetFinished();
        }
    }

    IteratorSimple& operator++() {
        this->SetFinished();
        return *this;
    }

    IteratorSimple operator++(int) {
        IteratorSimple iterator(*this);
        ++(*this);
        return iterator;
    }
};

}  // namespace improbable::phtree::v16

#endif  // PHTREE_V16_ITERATOR_SIMPLE_H
