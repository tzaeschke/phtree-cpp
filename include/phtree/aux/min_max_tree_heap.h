// Copyright Malte Skarupke 2020.
// Copyright Tilmann Zäschke 2023.
// Distributed under the Boost Software License, Version 1.0.
// (See http://www.boost.org/LICENSE_1_0.txt)

#ifndef PHTREE_AUX_MIN_MAX_TREE_HEAP_H
#define PHTREE_AUX_MIN_MAX_TREE_HEAP_H

#include "include/phtree/common/bpt_vector_tree.h"
#include "min_max_helpers.h"
#include <cstdint>
#include <functional>
#include <utility>

namespace phtree::aux {

namespace tree_heap {
template <typename It, typename Compare>
bool is_minmax_heap(It begin, It end, Compare&& compare) {
    uint64_t length = static_cast<uint64_t>(end - begin);
    auto test_index = [](uint64_t index, auto compare_index) {
        uint64_t first_child = minmax_heap_helpers::first_child_index(index);
        uint64_t second_child = first_child + 1;
        uint64_t first_grandchild = minmax_heap_helpers::first_child_index(first_child);
        uint64_t second_grandchild = first_grandchild + 1;
        uint64_t third_grandchild = minmax_heap_helpers::first_child_index(second_child);
        uint64_t fourth_grandchild = third_grandchild + 1;
        return compare_index(first_child) && compare_index(second_child) &&
            compare_index(first_grandchild) && compare_index(second_grandchild) &&
            compare_index(third_grandchild) && compare_index(fourth_grandchild);
    };
    for (uint64_t i = 0; i < length; ++i) {
        if (minmax_heap_helpers::is_min_item(i)) {
            auto compare_one = [&](uint64_t child) {
                return child >= length || !compare(begin[child], begin[i]);
            };
            if (!test_index(i, compare_one))
                return false;
        } else {
            auto compare_one = [&](uint64_t child) {
                return child >= length || !compare(begin[i], begin[child]);
            };
            if (!test_index(i, compare_one))
                return false;
        }
    }
    return true;
}

template <typename Data, typename Compare>
void push_minmax_heap(Data& data, Compare&& compare) {
    uint64_t length = static_cast<uint64_t>(data.size());
    uint64_t index = length - 1;
    uint64_t parent = minmax_heap_helpers::parent_index(index);
    auto value = std::move(data.back());
    if (minmax_heap_helpers::is_new_item_min(length)) {
        if (index == 0)
            static_cast<void>(0);
        else if (compare(data[parent], value)) {
            data[index] = std::move(data[parent]);
            index = parent;
            goto push_up_max;
        } else {
            for (;;) {
                {
                    uint64_t grandparent = minmax_heap_helpers::grandparent_index(index);
                    if (compare(value, data[grandparent])) {
                        data[index] = std::move(data[grandparent]);
                        index = grandparent;
                    } else
                        break;
                }
            push_up_min:
                if (!index)
                    break;
            }
        }
    } else if (compare(value, data[parent])) {
        data[index] = std::move(data[parent]);
        index = parent;
        goto push_up_min;
    } else {
    push_up_max:
        while (index > 2) {
            uint64_t grandparent = minmax_heap_helpers::grandparent_index(index);
            if (compare(data[grandparent], value)) {
                data[index] = std::move(data[grandparent]);
                index = grandparent;
            } else
                break;
        }
    }
    data[index] = std::move(value);
}

template <typename Data, typename Compare>
void pop_minmax_heap_min(Data& data, Compare&& compare) {
    uint64_t length = static_cast<uint64_t>(data.size()) - 1;
    if (length == 0)
        return;
    minmax_heap_helpers::push_down_min(
        data, std::exchange(data.back(), std::move(data.front())), 0, length, compare);
}

template <typename Data, typename Compare>
void pop_minmax_heap_max(Data& data, Compare&& compare) {
    uint64_t length = static_cast<uint64_t>(data.size()) - 1;
    if (length <= 1)
        return;

    uint64_t index = 1 + !!compare(data[1], data[2]);
    minmax_heap_helpers::push_down_max(
        data,
        std::exchange(data.back(), std::move(data[index])),
        index,
        length,
        std::forward<Compare>(compare));
}

template <typename It, typename Compare>
void make_minmax_heap(It begin, It end, Compare&& compare) {
    uint64_t length = end - begin;
    uint64_t index = length / 2;
    if (index == 0)
        return;
    // optimization: there can be only one item that has only one child
    // handling that item up front simplifies the second loop a little, since
    // we know that all other items have two children
    if ((length & 1) == 0) {
        --index;
        if (minmax_heap_helpers::is_min_item(index))
            minmax_heap_helpers::push_down_min_one_child_only(begin, index, compare);
        else
            minmax_heap_helpers::push_down_max_one_child_only(begin, index, compare);
        if (index == 0)
            return;
    }
    // optimization: half of all the items will have no grandchildren. this
    // simplifies the push_down function a lot, so we handle these items
    // first. we could then do another optimization where we know that
    // after the first half, the next quarter of items has grandchildren but
    // no great-grandchildren, but the code is already big enough
    if (length != 4) {
        uint64_t lowest_index_with_no_grandchildren = length / 4;
        for (;;) {
            int highest_bit = minmax_heap_helpers::highest_set_bit(index);
            uint64_t loop_until = std::max(
                lowest_index_with_no_grandchildren, (static_cast<uint64_t>(1) << highest_bit) - 1);
            --index;
            if (highest_bit & 1) {
                for (;; --index) {
                    minmax_heap_helpers::push_down_max_one_level_only(begin, index, compare);
                    if (index == loop_until)
                        break;
                }
            } else {
                for (;; --index) {
                    minmax_heap_helpers::push_down_min_one_level_only(begin, index, compare);
                    if (index == loop_until)
                        break;
                }
                if (index == 0)
                    return;
            }
            if (index == lowest_index_with_no_grandchildren)
                break;
        }
    }
    int highest_bit = minmax_heap_helpers::highest_set_bit(index);
    uint64_t loop_until = (static_cast<uint64_t>(1) << highest_bit) - 1;
    switch (highest_bit & 1) {
        for (;;) {
        case 0:
            for (;;) {
                --index;
                minmax_heap_helpers::push_down_min(
                    begin, std::move(begin[index]), index, length, compare);
                if (index == loop_until)
                    break;
            }
            if (index == 0)
                return;
            loop_until /= 2;
            [[fallthrough]];
        case 1:
            for (;;) {
                --index;
                minmax_heap_helpers::push_down_max(
                    begin, std::move(begin[index]), index, length, compare);
                if (index == loop_until)
                    break;
            }
            loop_until /= 2;
        }
    }
}
}  // namespace tree_heap

/**
 * A min-max heap that uses a vector-tree as underlying data structure.
 */
template <typename T, typename Compare>
class min_max_tree_heap {
    // TODO
    //   - reserve()
    //   - remove size() ?

    struct SwapComp {
        Compare comp;
        // template<class T>
        bool operator()(T const& x, T const& y) const {
            return !comp(x, y);
        }
    };

  public:
    const T& top() const {
        return arr[0];
    }

    const T& top_max() const {
        assert(arr.size() >= 3);
        uint64_t index = 1 + !!cmp(arr[1], arr[2]);
        return arr[index];
    }

    T& top_max() {
        assert(arr.size() >= 3);
        uint64_t index = 1 + !!cmp(arr[1], arr[2]);
        return arr[index];
    }

    template <typename... Args>
    void emplace(Args&&... args) {
        reserve();
        arr.emplace_back(std::forward<Args>(args)...);
        tree_heap::push_minmax_heap(arr, cmp);
    }

    void emplace(T&& x) {
        reserve();
        arr.emplace_back(std::move(x));
        tree_heap::push_minmax_heap(arr, cmp);
    }

    void emplace(const T& x) {
        reserve();
        arr.emplace_back(x);
        tree_heap::push_minmax_heap(arr, cmp);
    }

    void pop() {
        assert(!arr.empty());
        tree_heap::pop_minmax_heap_min(arr, cmp);
        arr.erase_back();
    }

    void pop_max() {
        assert(!arr.empty());
        tree_heap::pop_minmax_heap_max(arr, cmp);
        arr.erase_back();
    }

    [[nodiscard]] bool empty() const {
        return arr.empty();
    }

    [[nodiscard]] size_t size() const {
        return arr.size();
    }

  private:
    void reserve() {
        //        if (arr.capacity() == arr.size()) {
        //            arr.reserve(arr.capacity() * 2);
        //        }
    }
    // std::vector<T> arr;  // The heap array.
    ::phtree::bptree::detail::vector_tree<T, 64> arr;  // The heap array.
    // Compare cmp{};
    SwapComp cmp{};  // TODO create on the fly only?
};
}  // namespace phtree::aux

#endif  // PHTREE_AUX_MIN_MAX_TREE_HEAP_H
