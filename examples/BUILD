load("@rules_cc//cc:defs.bzl", "cc_binary")

package(default_visibility = ["//visibility:private"])

cc_binary(
    name = "example",
    srcs = ["example.cc"],
    visibility = [
        "//visibility:public",
    ],
    deps = [
        "//:phtree",
    ],
)
