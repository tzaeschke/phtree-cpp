cmake_minimum_required(VERSION 3.14)

# set the project name
project(phtree VERSION 1.3.0
        DESCRIPTION "PH-Tree C++"
        LANGUAGES CXX)


cmake_policy(SET CMP0077 NEW)

# ---------------------------------------------------------------------------------------
# Set default build to release
# ---------------------------------------------------------------------------------------
if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose Release or Debug" FORCE)
endif ()


# ---------------------------------------------------------------------------------------
# Build options
# ---------------------------------------------------------------------------------------
# example options
option(PHTREE_BUILD_ALL "Build examples, tests and benchmarks" OFF)

# example options
option(PHTREE_BUILD_EXAMPLES "Build examples" OFF)
#option(PHTREE_BUILD_EXAMPLE_HO "Build header only example" OFF)

# testing options
option(PHTREE_BUILD_TESTS "Build tests" OFF)
#option(PHTREE_BUILD_TESTS_HO "Build tests using the header only version" OFF)
option(PHTREE_CODE_COVERAGE "Collect coverage from test library" OFF)
if (PHTREE_CODE_COVERAGE)
    set(PHTREE_BUILD_TESTS ON)
endif ()

# bench options
option(PHTREE_BUILD_BENCHMARKS "Build benchmarks (Requires https://github.com/google/benchmark.git to be installed)" OFF)


# ---------------------------------------------------------------------------------------
# Compiler config
# ---------------------------------------------------------------------------------------
find_program(CCACHE_FOUND ccache)
if (CCACHE_FOUND)
    message("CCACHE is found")
    set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE ccache)
    set_property(GLOBAL PROPERTY RULE_LAUNCH_LINK ccache)
else (CCACHE_FOUND)
    message("CCACHE is NOT found")
endif (CCACHE_FOUND)

# specify the C++ standard
if (NOT CMAKE_CXX_STANDARD)
    set(CMAKE_CXX_STANDARD 17)
    set(CMAKE_CXX_STANDARD_REQUIRED ON)
endif ()

if (MSVC)
    #set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /std:c++17 /Wall")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /std:c++17")
    set(CMAKE_CXX_FLAGS_RELEASE "/O2")

    if (PHTREE_BUILD_TESTS OR PHTREE_BUILD_ALL)
        add_compile_options(/bigobj)
    endif ()

    # For google benchmark
    if (PHTREE_BUILD_BENCHMARKS) # OR PHTREE_BUILD_ALL)
        # This still doesn't work. This also breaks gtest
        # See for example
        # https://stackoverflow.com/questions/55376111/how-to-build-and-link-google-benchmark-using-cmake-in-windows
        # https://github.com/google/benchmark/issues/1348
        # https://github.com/google/benchmark/issues/639
        # set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS TRUE)
        # set(BUILD_SHARED_LIBS TRUE) #=TRUE
        # set(BENCHMARK_DOWNLOAD_DEPENDENCIES on)
        # set(BENCHMARK_ENABLE_GTEST_TESTS OFF)
    endif ()
else ()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 -Wall -Werror")
    if (PHTREE_BUILD_BENCHMARKS)
        set(CMAKE_CXX_FLAGS_RELEASE "-O3 -mavx -pthread")
    else ()
        set(CMAKE_CXX_FLAGS_RELEASE "-O3 -mavx")
    endif ()
    if (PHTREE_CODE_COVERAGE)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 -Wall -Werror -Wa,-mbig-obj")
    endif ()
endif ()

# ---------------------------------------------------------------------------------------
# Build binaries
# ---------------------------------------------------------------------------------------
add_subdirectory(phtree)

if (PHTREE_BUILD_EXAMPLES OR PHTREE_BUILD_ALL)
    message(STATUS "Generating examples")
    add_subdirectory(examples)
endif ()

if ((PHTREE_BUILD_BENCHMARKS OR PHTREE_BUILD_ALL) AND NOT MSVC)
    message(STATUS "Generating benchmarks")
    add_subdirectory(benchmark)
endif ()

if (PHTREE_BUILD_TESTS OR PHTREE_BUILD_ALL)
    message(STATUS "Generating tests")
    if (FALSE)
        add_compile_definitions(GTEST_HAS_ABSL=0)
        add_compile_definitions(GTEST_OS_WINDOWS_MOBILE=0)
        if (MSVC)
            add_compile_definitions(GTEST_OS_WINDOWS_MINGW=0)
        endif ()
        add_compile_definitions(GTEST_OS_LINUX_ANDROID=0)
        if (LINUX)
            add_compile_definitions(GTEST_OS_LINUX=1)
        else ()
            add_compile_definitions(GTEST_OS_LINUX=0)
        endif ()
        add_compile_definitions(
                GTEST_OS_WINDOWS_MOBILE=0
                GTEST_OS_WINDOWS_PHONE=0
                GTEST_OS_WINDOWS_RT=0
                GTEST_OS_ESP8266=0
                GTEST_OS_XTENSA=0)
    endif ()

    enable_testing()
    include(GoogleTest)
    add_subdirectory(test)
endif ()
