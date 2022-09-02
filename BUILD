package(default_visibility = ["//visibility:public"])

licenses(["notice"])  # Apache 2.0

# Expose license for external usage through bazel.
exports_files([
    "LICENSE",
])

# Platform configuration definitions for select()

config_setting(
    name = "linux",
    constraint_values = ["@bazel_tools//platforms:linux"],
)

config_setting(
    name = "macos",
    constraint_values = ["@bazel_tools//platforms:osx"],
)

config_setting(
    name = "macos_not_ios",
    constraint_values = ["@bazel_tools//platforms:osx"],
)

config_setting(
    name = "windows",
    constraint_values = ["@bazel_tools//platforms:windows"],
)

config_setting(
    name = "windows_debug",
    constraint_values = ["@bazel_tools//platforms:windows"],
    values = {
        "compilation_mode": "dbg",
    },
)

config_setting(
    name = "windows_release",
    constraint_values = ["@bazel_tools//platforms:windows"],
    values = {
        "compilation_mode": "opt",
    },
)

config_setting(
    name = "windows-x86_64",
    constraint_values = ["@bazel_tools//platforms:windows"],
)

# Buildifier

sh_binary(
    name = "buildifier",
    srcs = select(
        {
            ":linux": ["@buildifier_linux//file"],
            ":macos": ["@buildifier_macos//file"],
            ":windows": ["@buildifier_windows//file"],
        },
    ),
)

# Aspect-based clang-format

filegroup(
    name = "dot_clang_format",
    srcs = [".clang-format"],
)
