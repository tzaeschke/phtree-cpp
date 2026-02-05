load("@rules_shell//shell:sh_binary.bzl", "sh_binary")

package(default_visibility = ["//visibility:private"])

sh_binary(
    name = "tsan",
    srcs = ["tsan.sh"],
    data = [
        "tsan-suppressions.txt",
    ],
)
