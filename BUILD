package(default_visibility = ["//visibility:public"])

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

#exports_files(["phtree"])

filegroup(
    name = "phtree",
    srcs = glob([  "phtree/*.h",
            "phtree/common/*.h",
            ]),
#    hdrs = glob([  "phtree/*.h",
#            "phtree/common/*.h",
#            ]),
#    srcs = ["phtfoo.cc"],
#    data = [
#        "//my_package:exported_testdata",
#        "//my_package:mygroup",
#    ],
)

#cc_library(
#    name = "gtest",
#    srcs = glob(
#        include = [
#            "googletest/src/*.cc",
#            "googletest/src/*.h",
#            "googletest/include/gtest/**/*.h",
#            "googlemock/src/*.cc",
#            "googlemock/include/gmock/**/*.h",
#        ],
#        exclude = [
#            "googletest/src/gtest-all.cc",
#            "googletest/src/gtest_main.cc",
#            "googlemock/src/gmock-all.cc",
#            "googlemock/src/gmock_main.cc",
#        ],
#    ),
#    hdrs = glob([
#        "googletest/include/gtest/*.h",
#        "googlemock/include/gmock/*.h",
#    ]),
#    copts = select({
#        ":qnx": [],
#        ":windows": [],
#        "//conditions:default": ["-mavx"],
#    }),
#    features = select({
#        ":windows": ["windows_export_all_symbols"],
#        "//conditions:default": [],
#    }),
#    includes = [
#        "googlemock",
#        "googlemock/include",
#        "googletest",
#        "googletest/include",
#    ],
#)