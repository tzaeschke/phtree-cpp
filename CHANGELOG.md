# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](http://semver.org/spec/v2.0.0.html).

## [Unreleased]
### Added
- Added custom bucket implementation (similar to std::unordered_set). This improves update performance by 5%-20%.
  [#44](https://github.com/tzaeschke/phtree-cpp/issues/44)
- Added `PhTree.relocate(old_key, new_key)` function. This is **a lot faster** than using other methods.
  Added `PhTreeMultimap.relocate_all(old_key, new_key)` to move multiple values at once.
  [#43](https://github.com/tzaeschke/phtree-cpp/issues/43)
- Added try_emplace(key, value) and try_emplace(iter_hint, key, value)
  [#40](https://github.com/tzaeschke/phtree-cpp/issues/40)
- Added FilterBoxAABB and FilterSphereAABB as examples for filtering a PH-Tree with box keys
  [#33](https://github.com/tzaeschke/phtree-cpp/issues/33)
### Changed
- Simplified internals of erase(). [#47](https://github.com/tzaeschke/phtree-cpp/pull/47)
- Removed internal use of `std::optional()` to slightly reduce memory overhead
  [#38](https://github.com/tzaeschke/phtree-cpp/issues/38)
- Removed restrictions on bazel version [#35](https://github.com/tzaeschke/phtree-cpp/issues/35)
- **API BREAKING CHANGE**: API of filters have been changed to be more correct, explicit and flexible.
  [#21](https://github.com/tzaeschke/phtree-cpp/issues/21)
  - Correctness: Converters and distance functions are not copied unnecessarily anymore.
  - Explicit: 
    Filters *must* have a mandatory parameter for a converter reference. This ensures that the correct
            converter is used, probably `tree.converter()`.
  - Flexible: 
    Distance functions can be provided through a universal reference (forwarding reference).
    Also, filters are now movable and copyable. 
  
- **API BREAKING CHANGE**: Allow filtering on buckets in multimaps. Multimap filters have different functions
  and function signatures than normal `PhTree` filters. [#26](https://github.com/tzaeschke/phtree-cpp/issues/26)

## [1.2.0] - 2022-04-14
### Changed
- Bugfix: FilterSphere was not working correctly. [#27](https://github.com/tzaeschke/phtree-cpp/issues/27)
- Potentially **BREAKING CHANGE**: Refactored API of all methods that accept callbacks and filters to
  accept universal/forwarding references.
  Also changed filters and callback to not require `const` methods. 
  [#22](https://github.com/tzaeschke/phtree-cpp/issues/22)
- Clean up iterator implementations. [#19](https://github.com/tzaeschke/phtree-cpp/issues/19)
- Make PhTree and PhTreeMultimap movable (move-assign/copy). [#18](https://github.com/tzaeschke/phtree-cpp/issues/18)
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


[Unreleased]: https://github.com/improbable-eng/phtree-cpp/compare/v1.2.0...HEAD
[1.2.0]: https://github.com/improbable-eng/phtree-cpp/compare/v1.2.0...v1.1.0
[1.1.1]: https://github.com/improbable-eng/phtree-cpp/compare/v1.1.0...v1.1.1
[1.1.0]: https://github.com/improbable-eng/phtree-cpp/compare/v1.0.0...v1.1.0
[1.0.1]: https://github.com/improbable-eng/phtree-cpp/compare/v1.0.0...v1.0.1
[1.0.0]: https://github.com/improbable-eng/phtree-cpp/compare/v0.1.0...v1.0.0
