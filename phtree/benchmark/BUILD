package(default_visibility = ["//visibility:private"])

cc_library(
    name = "benchmark",
    testonly = True,
    srcs = [
        "logging.cc",
    ],
    hdrs = [
        "benchmark_util.h",
        "logging.h",
    ],
    visibility = [
        "//visibility:public",
    ],
    deps = [
        "@gbenchmark//:benchmark",
        "@spdlog",
    ],
    alwayslink = 1,
)

cc_binary(
    name = "count_mm_d_benchmark",
    testonly = True,
    srcs = [
        "count_mm_d_benchmark.cc",
    ],
    linkstatic = True,
    deps = [
        "//phtree",
        "//phtree/benchmark",
        "@gbenchmark//:benchmark",
        "@spdlog",
    ],
)

cc_binary(
    name = "erase_benchmark",
    testonly = True,
    srcs = [
        "erase_benchmark.cc",
    ],
    linkstatic = True,
    deps = [
        "//phtree",
        "//phtree/benchmark",
        "@gbenchmark//:benchmark",
        "@spdlog",
    ],
)

cc_binary(
    name = "erase_d_benchmark",
    testonly = True,
    srcs = [
        "erase_d_benchmark.cc",
    ],
    linkstatic = True,
    deps = [
        "//phtree",
        "//phtree/benchmark",
        "@gbenchmark//:benchmark",
        "@spdlog",
    ],
)

cc_binary(
    name = "extent_benchmark",
    testonly = True,
    srcs = [
        "extent_benchmark.cc",
    ],
    linkstatic = True,
    deps = [
        "//phtree",
        "//phtree/benchmark",
        "@gbenchmark//:benchmark",
        "@spdlog",
    ],
)

cc_binary(
    name = "extent_benchmark_weird",
    testonly = True,
    srcs = [
        "extent_benchmark_weird.cc",
    ],
    linkstatic = True,
    deps = [
        "//phtree",
        "//phtree/benchmark",
        "@gbenchmark//:benchmark",
        "@spdlog",
    ],
)

cc_binary(
    name = "find_benchmark",
    testonly = True,
    srcs = [
        "find_benchmark.cc",
    ],
    linkstatic = True,
    deps = [
        "//phtree",
        "//phtree/benchmark",
        "@gbenchmark//:benchmark",
        "@spdlog",
    ],
)

cc_binary(
    name = "insert_benchmark",
    testonly = True,
    srcs = [
        "insert_benchmark.cc",
    ],
    linkstatic = True,
    deps = [
        "//phtree",
        "//phtree/benchmark",
        "@gbenchmark//:benchmark",
        "@spdlog",
    ],
)

cc_binary(
    name = "insert_d_benchmark",
    testonly = True,
    srcs = [
        "insert_d_benchmark.cc",
    ],
    linkstatic = True,
    deps = [
        "//phtree",
        "//phtree/benchmark",
        "@gbenchmark//:benchmark",
        "@spdlog",
    ],
)

cc_binary(
    name = "insert_box_d_benchmark",
    testonly = True,
    srcs = [
        "insert_box_d_benchmark.cc",
    ],
    linkstatic = True,
    deps = [
        "//phtree",
        "//phtree/benchmark",
        "@gbenchmark//:benchmark",
        "@spdlog",
    ],
)

cc_binary(
    name = "knn_d_benchmark",
    testonly = True,
    srcs = [
        "knn_d_benchmark.cc",
    ],
    linkstatic = True,
    deps = [
        "//phtree",
        "//phtree/benchmark",
        "@gbenchmark//:benchmark",
        "@spdlog",
    ],
)

cc_binary(
    name = "query_benchmark",
    testonly = True,
    srcs = [
        "query_benchmark.cc",
    ],
    linkstatic = True,
    deps = [
        "//phtree",
        "//phtree/benchmark",
        "@gbenchmark//:benchmark",
        "@spdlog",
    ],
)

cc_binary(
    name = "query_box_d_benchmark",
    testonly = True,
    srcs = [
        "query_box_d_benchmark.cc",
    ],
    linkstatic = True,
    deps = [
        "//phtree",
        "//phtree/benchmark",
        "@gbenchmark//:benchmark",
        "@spdlog",
    ],
)

cc_binary(
    name = "query_d_benchmark",
    testonly = True,
    srcs = [
        "query_d_benchmark.cc",
    ],
    linkstatic = True,
    deps = [
        "//phtree",
        "//phtree/benchmark",
        "@gbenchmark//:benchmark",
        "@spdlog",
    ],
)

cc_binary(
    name = "query_mm_d_benchmark",
    testonly = True,
    srcs = [
        "query_mm_d_benchmark.cc",
    ],
    linkstatic = True,
    deps = [
        "//phtree",
        "//phtree/benchmark",
        "@gbenchmark//:benchmark",
        "@spdlog",
    ],
)

cc_binary(
    name = "query_mm_box_d_benchmark",
    testonly = True,
    srcs = [
        "query_mm_box_d_benchmark.cc",
    ],
    linkstatic = True,
    deps = [
        "//phtree",
        "//phtree/benchmark",
        "@gbenchmark//:benchmark",
        "@spdlog",
    ],
)

cc_binary(
    name = "update_d_benchmark",
    testonly = True,
    srcs = [
        "update_d_benchmark.cc",
    ],
    linkstatic = True,
    deps = [
        "//phtree",
        "//phtree/benchmark",
        "@gbenchmark//:benchmark",
        "@spdlog",
    ],
)

cc_binary(
    name = "update_mm_d_benchmark",
    testonly = True,
    srcs = [
        "update_mm_d_benchmark.cc",
    ],
    linkstatic = True,
    deps = [
        "//phtree",
        "//phtree/benchmark",
        "@gbenchmark//:benchmark",
        "@spdlog",
    ],
)

cc_binary(
    name = "update_mm_box_d_benchmark",
    testonly = True,
    srcs = [
        "update_mm_box_d_benchmark.cc",
    ],
    linkstatic = True,
    deps = [
        "//phtree",
        "//phtree/benchmark",
        "@gbenchmark//:benchmark",
        "@spdlog",
    ],
)

cc_binary(
    name = "update_box_d_benchmark",
    testonly = True,
    srcs = [
        "update_box_d_benchmark.cc",
    ],
    linkstatic = True,
    deps = [
        "//phtree",
        "//phtree/benchmark",
        "@gbenchmark//:benchmark",
        "@spdlog",
    ],
)
