package(default_visibility = ["//visibility:public"])

cc_library(
    name = "phtree",
    srcs = glob(
        include = [
            "common/*.h",
            "v16/*.h",
        ],
    ),
    hdrs = [
        "phtree.h",
        "phtree_multimap.h",
    ],
    linkstatic = True,
    visibility = [
        "//visibility:public",
    ],
    deps = [
        "//phtree/v16",
    ],
)
