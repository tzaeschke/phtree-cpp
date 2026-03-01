
# Create Release (document WIP)

## Verify code

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

```bash
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

#### Start

Checkout personal fork of bazel-central-registry.
Create branch `phtree-cpp@1.7.0`.

#### After modifying any file

Recreate hashes for `source.json`
```shell
bazel run -- //tools:update_integrity phtree-cpp
```
 This also accepts `--version` but defaults to latest version.

#### Validate 

```shell
bazel run -- //tools:bcr_validation --check=phtree-cpp@1.7.0
```

#### Test locally

```bash
bazel run //tools:setup_presubmit_repos -- --module phtree-cpp@1.7.0
```

There is also a GitHub action for release automation:
https://github.com/bazel-contrib/publish-to-bcr


#### Test with Test repos

Check out `test-phtree-cpp-bazel` and run:

```shell
bazel shutdown && bazel build --enable_bzlmod --registry="file:///home/.../bazel-central-registry" --lockfile_mode=off @phtree-cpp//:phtree-cpp
--> Target might be incorrect
```