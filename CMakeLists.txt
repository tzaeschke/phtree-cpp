cmake_minimum_required(VERSION 3.14)

# set the project name
project(PH_Tree_Main VERSION 1.2.0
        DESCRIPTION "PH-Tree C++"
        LANGUAGES CXX)

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
option(PHTREE_BUILD_EXAMPLES "Build examples" OFF)
#option(PHTREE_BUILD_EXAMPLE_HO "Build header only example" OFF)

# testing options
option(PHTREE_BUILD_TESTS "Build tests" OFF)
#option(PHTREE_BUILD_TESTS_HO "Build tests using the header only version" OFF)

# bench options
option(PHTREE_BUILD_BENCHMARKS "Build benchmarks (Requires https://github.com/google/benchmark.git to be installed)" OFF)

# ---------------------------------------------------------------------------------------
# Compiler config
# ---------------------------------------------------------------------------------------
# specify the C++ standard
if (NOT CMAKE_CXX_STANDARD)
    set(CMAKE_CXX_STANDARD 17)
    set(CMAKE_CXX_STANDARD_REQUIRED ON)
endif ()

if (MSVC)
    #set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /std:c++17 /Wall")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /std:c++17")
    set(CMAKE_CXX_FLAGS_RELEASE "/O2")

    ##    set(CMAKE_CXX_FLAGS "-DNOMINMAX ${CMAKE_CXX_FLAGS}") # exclude M$ min/max macros
    ##    set(CMAKE_CXX_FLAGS "/wd4996 ${CMAKE_CXX_FLAGS}") # don't warn about use of plain C functions without (non-portable) "_s" suffix
    ##    set(CMAKE_EXE_LINKER_FLAGS "/WX:NO ${CMAKE_EXE_LINKER_FLAGS}" ) # don't treat warnings as compile errors--gtest doesn't build
    ##    #set(CMAKE_CXX_FLAGS_DEBUG "/analyze ${CMAKE_CXX_FLAGS_DEBUG}")
else ()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 -Wall -Werror")
    if (PHTREE_BUILD_BENCHMARKS)
        set(CMAKE_CXX_FLAGS_RELEASE "-O3 -mavx -pthread")
    else ()
        set(CMAKE_CXX_FLAGS_RELEASE "-O3 -mavx")
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

if (PHTREE_BUILD_BENCHMARKS OR PHTREE_BUILD_ALL)
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
