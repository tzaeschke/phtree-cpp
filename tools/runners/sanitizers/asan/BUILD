load("@rules_shell//shell:sh_binary.bzl", "sh_binary")

package(default_visibility = ["//visibility:private"])

sh_binary(
    name = "asan",
    srcs = ["asan.sh"],
    data = [
        "asan-suppressions.txt",
        "lsan-suppressions.txt",
    ],
)
