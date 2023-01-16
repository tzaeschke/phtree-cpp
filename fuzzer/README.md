Fuzzing
=======

Requirements:
 * `clang`.
 * libFuzzer: https://github.com/google/fuzzing/blob/master/tutorial/libFuzzerTutorial.md 

```
CC=clang bazel run //fuzzer:b_plus_multimap_fuzzer --config=fuzz
```