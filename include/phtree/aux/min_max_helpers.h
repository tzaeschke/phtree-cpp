// Copyright Malte Skarupke 2020.
// Distributed under the Boost Software License, Version 1.0.
// (See http://www.boost.org/LICENSE_1_0.txt)

#ifndef PHTREE_AUX_MIN_MAX_HELPERS_H
#define PHTREE_AUX_MIN_MAX_HELPERS_H

namespace phtree::aux::minmax_heap_helpers {
// returns the index of the highest set bit. undefined if no bits are set.
// examples:
// highest_set_bit(1) = 0
// highest_set_bit(4) = 2
// highest_set_bit(55) = 5
inline int highest_set_bit(uint64_t i) {
#ifdef _MSC_VER
    unsigned long result;
    _BitScanReverse64(&result, i);
    return result;
#else
    return 63 - __builtin_clzl(i);
#endif
}

inline bool is_new_item_min(uint64_t length) {
    return (highest_set_bit(length) & 1) == 0;
}

inline bool is_min_item(uint64_t index) {
    return is_new_item_min(index + 1);
}

inline uint64_t grandparent_index(uint64_t index) {
    return (index - 3) / 4;
}

inline uint64_t parent_index(uint64_t index) {
    return (index - 1) / 2;
}

inline uint64_t first_child_index(uint64_t index) {
    return (index * 2) + 1;
}
inline uint64_t last_grandchild_index(uint64_t index) {
    return (index * 4) + 6;
}
template <typename Data, typename Compare>
uint64_t smallest_descendant(
    Data& data,
    uint64_t length,
    uint64_t first_child,
    uint64_t first_grandchild,
    Compare&& compare) {
    uint64_t second_child = first_child + 1;
    if (first_grandchild >= length)
        return first_child +
            (second_child != length && compare(data[second_child], data[first_child]));
    uint64_t second_grandchild = first_grandchild + 1;
    if (second_grandchild == length)
        return compare(data[first_grandchild], data[second_child]) ? first_grandchild
                                                                   : second_child;
    uint64_t min_grandchild =
        first_grandchild + !!compare(data[second_grandchild], data[first_grandchild]);
    uint64_t third_grandchild = second_grandchild + 1;
    if (third_grandchild == length)
        return compare(data[min_grandchild], data[second_child]) ? min_grandchild : second_child;
    else
        return compare(data[min_grandchild], data[third_grandchild]) ? min_grandchild
                                                                     : third_grandchild;
}
template <typename Data, typename Compare>
uint64_t largest_descendant(
    Data& data,
    uint64_t length,
    uint64_t first_child,
    uint64_t first_grandchild,
    Compare&& compare) {
    uint64_t second_child = first_child + 1;
    if (first_grandchild >= length)
        return first_child +
            (second_child != length && compare(data[first_child], data[second_child]));
    uint64_t second_grandchild = first_grandchild + 1;
    if (second_grandchild == length)
        return compare(data[second_child], data[first_grandchild]) ? first_grandchild
                                                                   : second_child;
    uint64_t max_grandchild =
        first_grandchild + !!compare(data[first_grandchild], data[second_grandchild]);
    uint64_t third_grandchild = second_grandchild + 1;
    if (third_grandchild == length)
        return compare(data[second_child], data[max_grandchild]) ? max_grandchild : second_child;
    else
        return compare(data[max_grandchild], data[third_grandchild]) ? third_grandchild
                                                                     : max_grandchild;
}

template <typename Data, typename Compare>
void push_down_min(
    Data& data,
    typename Data::value_type value,
    uint64_t index,
    uint64_t length,
    Compare&& compare) {
    using std::swap;
    for (;;) {
        uint64_t last_grandchild = last_grandchild_index(index);
        if (last_grandchild < length) {
            // auto it = data.begin() + last_grandchild; // TODO?
            auto it = last_grandchild;
            uint64_t min_first_half = last_grandchild - 2 - !!compare(data[it - 3], data[it - 2]);
            uint64_t min_second_half = last_grandchild - !!compare(data[it - 1], data[it]);
            uint64_t smallest = compare(data[min_second_half], data[min_first_half])
                ? min_second_half
                : min_first_half;
            if (!compare(data[smallest], value))
                break;
            data[index] = std::move(data[smallest]);
            index = smallest;
            uint64_t parent = parent_index(index);
            if (compare(data[parent], value))
                swap(data[parent], value);
        } else {
            uint64_t first_child = first_child_index(index);
            if (first_child >= length)
                break;
            uint64_t first_grandchild = last_grandchild - 3;
            uint64_t smallest =
                smallest_descendant(data, length, first_child, first_grandchild, compare);
            if (!compare(data[smallest], value))
                break;
            data[index] = std::move(data[smallest]);
            index = smallest;
            if (smallest < first_grandchild)
                break;
            uint64_t parent = parent_index(index);
            if (compare(data[parent], value)) {
                data[index] = std::move(data[parent]);
                index = parent;
            }
            break;
        }
    }
    data[index] = std::move(value);
}

template <typename It, typename Compare>
void push_down_min_one_child_only(It begin, uint64_t index, Compare&& compare) {
    using std::swap;
    uint64_t child = first_child_index(index);
    if (compare(begin[child], begin[index]))
        swap(begin[index], begin[child]);
}

template <typename It, typename Compare>
void push_down_min_one_level_only(It begin, uint64_t index, Compare&& compare) {
    using std::swap;
    uint64_t first_child = first_child_index(index);
    uint64_t smaller_child = first_child + !!compare(begin[first_child + 1], begin[first_child]);
    if (compare(begin[smaller_child], begin[index]))
        swap(begin[index], begin[smaller_child]);
}

template <typename Data, typename Compare>
void push_down_max(
    Data& data,
    typename Data::value_type value,
    uint64_t index,
    uint64_t length,
    Compare&& compare) {
    using std::swap;
    for (;;) {
        uint64_t last_grandchild = last_grandchild_index(index);
        if (last_grandchild < length) {
            // auto it = data.begin() + last_grandchild; // TODO?
            auto it = last_grandchild;
            uint64_t max_first_half = last_grandchild - 2 - !!compare(data[it - 2], data[it - 3]);
            uint64_t max_second_half = last_grandchild - !!compare(data[it], data[it - 1]);
            uint64_t largest = compare(data[max_first_half], data[max_second_half])
                ? max_second_half
                : max_first_half;
            if (!compare(value, data[largest]))
                break;
            data[index] = std::move(data[largest]);
            index = largest;
            uint64_t parent = parent_index(index);
            if (compare(value, data[parent]))
                swap(data[parent], value);
        } else {
            uint64_t first_child = first_child_index(index);
            if (first_child >= length)
                break;
            uint64_t first_grandchild = last_grandchild - 3;
            uint64_t largest =
                largest_descendant(data, length, first_child, first_grandchild, compare);
            if (!compare(value, data[largest]))
                break;
            data[index] = std::move(data[largest]);
            index = largest;
            if (largest < first_grandchild)
                break;
            uint64_t parent = parent_index(index);
            if (compare(value, data[parent])) {
                data[index] = std::move(data[parent]);
                index = parent;
            }
            break;
        }
    }
    data[index] = std::move(value);
}

template <typename It, typename Compare>
void push_down_max_one_child_only(It begin, uint64_t index, Compare&& compare) {
    using std::swap;
    uint64_t child = first_child_index(index);
    if (compare(begin[index], begin[child]))
        swap(begin[index], begin[child]);
}

template <typename It, typename Compare>
void push_down_max_one_level_only(It begin, uint64_t index, Compare&& compare) {
    using std::swap;
    uint64_t first_child = first_child_index(index);
    uint64_t bigger_child = first_child + !!compare(begin[first_child], begin[first_child + 1]);
    if (compare(begin[index], begin[bigger_child]))
        swap(begin[index], begin[bigger_child]);
}

}  // namespace minmax_heap_helpers


#endif  // PHTREE_AUX_MIN_MAX_HELPERS_H
