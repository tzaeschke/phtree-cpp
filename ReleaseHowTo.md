



https://github.com/tzaeschke/phtree-cpp/releases/v1.6.2/download/v1.6.2.tar.gz


## Create Release


### GitHub

Create release.


### Archive

```bash
git archive --format=zip v1.6.2 > phtree-cpp-v1.6.2.zip
```

Upload to GitHub release.

Resulting link:
```
https://github.com/tzaeschke/phtree-cpp/releases/download/v1.6.2/phtree-cpp-v1.6.2.zip
```


### Bazel Central Repository

Documentation: 
https://github.com/bazelbuild/bazel-central-registry/blob/main/docs/README.md


#### Start

Checkout personal fork of bazel-central-registry.
Create branch phtree-cpp@1.6.2

#### After modifying any file

Recreate hashes for `source.json`
```shell
bazel run -- //tools:update_integrity phtree-cpp
```
 This also accepts `--version` but defaults to latest version.

#### Validate 

```shell
bazel run -- //tools:bcr_validation --check=phtree-cpp@1.6.2
```

#### Test locally

```bash
bazel run //tools:setup_presubmit_repos -- --module phtree-cpp@1.6.2
```

There is also a GitHub action for release automation:
https://github.com/bazel-contrib/publish-to-bcr