load("@buildifier_prebuilt//:rules.bzl", "buildifier")
load("@rules_cc//cc:defs.bzl", "cc_library")

package(default_visibility = ["//visibility:public"])

licenses(["notice"])  # Apache 2.0

# Expose license for external usage through bazel.
exports_files([
    "LICENSE",
])

# Platform configuration definitions for select()

config_setting(
    name = "linux",
    constraint_values = ["@platforms//os:linux"],
)

config_setting(
    name = "macos",
    constraint_values = ["@platforms//os:osx"],
)

config_setting(
    name = "macos_not_ios",
    constraint_values = ["@platforms//os:osx"],
)

config_setting(
    name = "windows",
    constraint_values = ["@platforms//os:windows"],
)

config_setting(
    name = "windows_debug",
    constraint_values = ["@platforms//os:windows"],
    values = {
        "compilation_mode": "dbg",
    },
)

config_setting(
    name = "windows_release",
    constraint_values = ["@platforms//os:windows"],
    values = {
        "compilation_mode": "opt",
    },
)

config_setting(
    name = "windows-x86_64",
    constraint_values = ["@platforms//os:windows"],
)

# Buildifier
buildifier(
    name = "buildifier.fix",
    exclude_patterns = ["./.git/*"],
    lint_mode = "fix",
    mode = "fix",
)

buildifier(
    name = "buildifier.check",
    exclude_patterns = [
        "./.git/*",
    ],
    lint_mode = "warn",
    mode = "diff",
)

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
