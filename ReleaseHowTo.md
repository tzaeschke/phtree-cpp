
# Create Release (document WIP)

## Verify code

### Run tests
```shell
CC=gcc bazel test ... --config=asan
CC=clang bazel test ... --config=asan
CC=clang bazel test ... --config=msan
CC=clang bazel test ... --config=ubsan
```

Record versions: e.g. gcc 13.3, clang 18.1.3

```shell
bazel run //:buildifier.check
bazel run //:buildifier.fix
```

### Run examples
```shell
bazel run //examples:example
```

### Test release with
- test-phtree-cpp-bazel
- test-phtree-cpp-cmake

### Compilers

- Test with gcc + clang on Linux
- Test with VisualStudio on Windows
- Update references of tested compilers in README.md

## Release Prep

* Create branch/PR `release/1.7.0`
* Update version in 
  * MODULE.bazel
  * CMakeLists.txt
  * README.md
  * CHANGELOG.md (don't forget links at the bottom)
* Commit, push, merge

## Create GitHub release 

...

## Create Archive

Checkout main branch and pull the tag.

```bash
git checkout main
git pull
git archive --format=zip v1.7.0 > phtree-cpp-v1.7.0.zip
```

Upload to GitHub release.

Resulting link:
```
https://github.com/tzaeschke/phtree-cpp/releases/download/v1.7.0/phtree-cpp-v1.7.0.zip
```
or 
```
https://github.com/tzaeschke/phtree-cpp/releases/v1.7.0/download/v1.7.0.tar.gz
```


## Bazel Central Repository -- Unfinished documentation

Documentation: 
https://github.com/bazelbuild/bazel-central-registry/blob/main/docs/README.md


PH-Tree build target:
```
bazel build //:phtree
```

#### Bazel Central Repository Update

1. Checkout personal fork of bazel-central-registry.
   Create branch `phtree-cpp@1.7.0`.

2. **Initially only**: #### Initially
   Only for first upload: `bazel run //tools:add_module`
   Test target: `//...`


3. Patch MODULE.bazel:
   ```shell
   diff -ruN ../phtree-cpp/MODULE.bazel modules/phtree-cpp/1.7.0/MODULE.bazel > modules/phtree-cpp/1.7.0/patches/module_dot_bazel.patch 
   ```
   Update first two lines of patch file to
   ```
   --- MODULE.bazel
   +++ MODULE.bazel
   ```

4. Recreate hashes for `source.json`
   ```shell
   bazel run -- //tools:update_integrity phtree-cpp
   ```
   This also accepts `--version` but defaults to latest version.


5.  Validate 

   ```shell
   bazel run -- //tools:bcr_validation --check=phtree-cpp@1.7.0
   ```

6. Test locally

   1. Build local repo   
      ```bash
      bazel run //tools:setup_presubmit_repos -- --module phtree-cpp@1.7.0
      ```

      There is also a GitHub action for release automation:
      https://github.com/bazel-contrib/publish-to-bcr

   2. Test with Test repos

      Check out `test-phtree-cpp-bazel` and run (**replace path with actual path**):

      ```shell
      bazel shutdown && bazel build --enable_bzlmod --registry="file:///home/ FIX PATH /bazel-central-registry" --lockfile_mode=off ...
      ```
7. Commit and push
   ```
   git commit -m phtree-cpp@1.7.0
   git push 
   ```

### Original presubmit.yml file

The following has been removed for now because gtest causes compilation errors on windows.
```yaml
  tasks:
    run_test_module:
      name: Run test module
      platform: ${{ platform }}
      bazel: ${{ bazel }}
      build_targets:
      - '@phtree-cpp//test/...'
      test_targets:
      - '@phtree-cpp//test/...'
```
-> 
```shell
external/googletest~/googletest/src/gtest-filepath.cc(30): fatal error C1083: Cannot open include file: 'gtest/internal/gtest-filepath.h': No such file or directory
(15:48:44) ERROR: C:/b/he4gwbce/external/googletest~/BUILD.bazel:102:11: Compiling googlemock/src/gmock-matchers.cc failed: (Exit 2): cl.exe failed: error executing CppCompile command (from target @@googletest~//:gtest)
```