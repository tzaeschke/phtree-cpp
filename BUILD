package(default_visibility = ["//visibility:public"])

licenses(["notice"])  # Apache 2.0

# Expose license for external usage through bazel.
exports_files([
    "LICENSE",
])

# Aspect-based clang-format

filegroup(
    name = "dot_clang_format",
    srcs = [".clang-format"],
)

cc_library(
    name = "phtree",
    srcs = glob(
        include = [
            "include/**/*.h",
        ],
    ),
    hdrs = [
        "include/phtree/converter.h",
        "include/phtree/distance.h",
        "include/phtree/filter.h",
        "include/phtree/phtree.h",
        "include/phtree/phtree_multimap.h",
    ],
    includes = [
        "include",
    ],
    linkstatic = True,
    visibility = [
        "//visibility:public",
    ],
    deps = [
        "//include/phtree/common",
        "//include/phtree/v16",
    ],
)
