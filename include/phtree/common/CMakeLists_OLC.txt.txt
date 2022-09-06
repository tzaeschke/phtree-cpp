cmake_minimum_required(VERSION 3.14)

target_sources(phtree
        PRIVATE
        common.h
        base_types.h
        bits.h
        flat_array_map.h
        flat_sparse_map.h
        debug_helper.h
        tree_stats.h
        )
