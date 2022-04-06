cmake_minimum_required(VERSION 3.14)

target_sources(phtree
        PRIVATE
        debug_helper_v16.h
        node.h
        entry.h
        iterator_base.h
        iterator_full.h
        iterator_hc.h
        iterator_knn_hs.h
        iterator_with_parent.h
        phtree_v16.h
        )
