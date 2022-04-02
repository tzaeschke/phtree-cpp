# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](http://semver.org/spec/v2.0.0.html).

## [Unreleased]
### Changed
- Make PhTree and PhTreeMultimap moveable (move-assign/copy). [#18](https://github.com/tzaeschke/phtree-cpp/issues/18)
- Potentially **BREAKING CHANGE** when using `IsNodeValid()` in provided filters:
  Changed `bit_width_t` from `uin16_t` to `uint32_t`. This improves performance of 3D insert/emplace
  on small datasets by up to 15%. To avoid warnings that meant that the API of `FilterAABB` and `FilterSphere` 
  had to be changed to accept `uint32_t` instead of `int`. This may break some implementations.
  [#17](https://github.com/tzaeschke/phtree-cpp/pull/17)
- DIM>8 now uses custom b_plus_tree_map instead of std::map. This improves performance for all operations, e.g.
  window queries on large datasets are up to 4x faster. Benchmarks results can be found in the issue. 
  [#14](https://github.com/tzaeschke/phtree-cpp/issues/14)
- postfix/infix field moved from Node to Entry. This avoids indirections and improves performance of most by ~10%.
  operations by 5-15%.  [#11](https://github.com/tzaeschke/phtree-cpp/issues/11)
- Entries now use 'union' to store children.  [#9](https://github.com/tzaeschke/phtree-cpp/issues/9)
- Avoid unnecessary find() when removing a node. [#5](https://github.com/tzaeschke/phtree-cpp/issues/5)
- Avoid unnecessary key copy when inserting a node. [#4](https://github.com/tzaeschke/phtree-cpp/issues/4)
- for_each(callback, filter) was traversing too many nodes. [#2](https://github.com/tzaeschke/phtree-cpp/issues/2)
- Build improvements for bazel/cmake

## [1.1.1] - 2022-01-30
### Changed
- Replaced size() in filters with DIM [#26](https://github.com/improbable-eng/phtree-cpp/pull/26)

## [1.1.0] - 2022-01-25
### Added
- FilterSphere for filtering by sphere constraint (by ctbur) [#16](https://github.com/improbable-eng/phtree-cpp/pull/16)
- IEEE converter for 32bit float, see `distance.h` (by ctbur) [#18](https://github.com/improbable-eng/phtree-cpp/pull/18)

### Changed
- Performance improvement for updates and queries: removed use of `std::variant`. [#23](https://github.com/improbable-eng/phtree-cpp/pull/23)
- Fixed imports  `<climits>` -> `<limits>` (by ctbur) [#15](https://github.com/improbable-eng/phtree-cpp/pull/15)
- Cleaned up build scripts [#21](https://github.com/improbable-eng/phtree-cpp/pull/21)
- Fixed warnings: [#20](https://github.com/improbable-eng/phtree-cpp/pull/20)
  - "unused function argument" warnings
  - gcc/clang warnings
  - MSVC warnings
  - reserved identifier warnings (identifiers starting with `_`)
- typos in README.md [#22](https://github.com/improbable-eng/phtree-cpp/pull/22)

## [1.0.1] - 2021-05-06
### Changed
- replaced compilation flag `-fpermissive` with `-Werror`, and fixed all warnings/errors, see issue #10

## [1.0.0] - 2021-03-23
### Added
- API: `MultiMap`: A wrapper that makes PH-Tree behave as a multi-map.
- API: `erase(iterator)`
- API: `emplace_hint(iterator, ...)`
- API for `PhTreeF` and `PhTreeBoxF`: 32bit floating point options
- Support for custom key classes

### Changed
- BREAKING CHANGE: The query functions now require a query box as input (instead of a min/max point pair)
- BREAKING CHANGE: `phtree_box_d.h` has been removed, please use `phtree.h instead.
- BREAKING CHANGE: `phtree_d.h` has been removed, please use `phtree.h` instead.
- BREAKING CHANGE: Data converters (IEEE, Multiply, etc) are now structs i.o. functions/functors
- BREAKING CHANGE: `PhFilterNoOp` has been renamed to `FilterNoOp`
- BREAKING CHANGE: kNN queries now always require the distance function to be specified.
- BREAKING CHANGE: Preprocessors have been refactored and renamed to Converter/ScalarConverter
- Moved CI builds from Travis to GitHub actions

### Removed
- Nothing.

### Fixed
- GCC warnings from `-Wsign-compare` and `-Wsequence-point`.


## 0.1.0 - 2020-07-02
### Added
- Initial version.

### Changed
- Nothing.

### Removed
- Nothing.

### Fixed
- Nothing.


[Unreleased]: https://github.com/improbable-eng/phtree-cpp/compare/v1.1.1...HEAD
[1.1.1]: https://github.com/improbable-eng/phtree-cpp/compare/v1.1.0...v1.1.1
[1.1.0]: https://github.com/improbable-eng/phtree-cpp/compare/v1.0.0...v1.1.0
[1.0.1]: https://github.com/improbable-eng/phtree-cpp/compare/v1.0.0...v1.0.1
[1.0.0]: https://github.com/improbable-eng/phtree-cpp/compare/v0.1.0...v1.0.0
