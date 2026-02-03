load("@rules_shell//shell:sh_binary.bzl", "sh_binary")

package(default_visibility = ["//visibility:private"])

sh_binary(
    name = "msan",
    srcs = ["msan.sh"],
    data = [
        "msan-suppressions.txt",
    ],
)
