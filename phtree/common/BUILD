package(default_visibility = ["//visibility:private"])

cc_library(
    name = "common",
    hdrs = [
        "base_types.h",
        "bits.h",
        "common.h",
        "converter.h",
        "debug_helper.h",
        "distance.h",
        "filter.h",
        "flat_array_map.h",
        "b_plus_tree_map.h",
        "flat_sparse_map.h",
        "tree_stats.h",
    ],
    visibility = [
        "//visibility:public",
    ],
    deps = [
    ],
)

cc_test(
    name = "base_types_test",
    timeout = "long",
    srcs = [
        "base_types_test.cc",
    ],
    linkstatic = True,
    deps = [
        ":common",
        "//phtree/testing/gtest_main",
    ],
)

cc_test(
    name = "bits_test",
    timeout = "long",
    srcs = [
        "bits_test.cc",
    ],
    linkstatic = True,
    deps = [
        ":common",
        "//phtree/testing/gtest_main",
    ],
)

cc_test(
    name = "common_test",
    timeout = "long",
    srcs = [
        "common_test.cc",
    ],
    linkstatic = True,
    deps = [
        ":common",
        "//phtree/testing/gtest_main",
    ],
)

cc_test(
    name = "distance_test",
    timeout = "long",
    srcs = [
        "distance_test.cc",
    ],
    linkstatic = True,
    deps = [
        ":common",
        "//phtree/testing/gtest_main",
    ],
)

cc_test(
    name = "filter_test",
    timeout = "long",
    srcs = [
        "filter_test.cc",
    ],
    linkstatic = True,
    deps = [
        ":common",
        "//phtree/testing/gtest_main",
    ],
)

cc_test(
    name = "flat_array_map_test",
    timeout = "long",
    srcs = [
        "flat_array_map_test.cc",
    ],
    linkstatic = True,
    deps = [
        ":common",
        "//phtree/testing/gtest_main",
    ],
)

cc_test(
    name = "b_plus_tree_map_test",
    timeout = "long",
    srcs = [
        "b_plus_tree_map_test.cc",
    ],
    linkstatic = True,
    deps = [
        ":common",
        "//phtree/testing/gtest_main",
    ],
)

cc_test(
    name = "flat_sparse_map_test",
    timeout = "long",
    srcs = [
        "flat_sparse_map_test.cc",
    ],
    linkstatic = True,
    deps = [
        ":common",
        "//phtree/testing/gtest_main",
    ],
)

cc_test(
    name = "preprocessor_test",
    timeout = "long",
    srcs = [
        "converter_test.cc",
    ],
    linkstatic = True,
    deps = [
        ":common",
        "//phtree/testing/gtest_main",
    ],
)
