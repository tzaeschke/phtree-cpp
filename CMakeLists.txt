cmake_minimum_required(VERSION 3.14)

# set the project name
project(PH_Tree_Main VERSION 1.2.0
        DESCRIPTION "PH-Tree C++"
        LANGUAGES CXX)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# specify the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)
if(MSVC)
    #set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /std:c++17 /Wall")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /std:c++17")
    set(CMAKE_CXX_FLAGS_RELEASE "/O2")

##    set(CMAKE_CXX_FLAGS "-DNOMINMAX ${CMAKE_CXX_FLAGS}") # exclude M$ min/max macros
##    set(CMAKE_CXX_FLAGS "/wd4996 ${CMAKE_CXX_FLAGS}") # don't warn about use of plain C functions without (non-portable) "_s" suffix
##    set(CMAKE_EXE_LINKER_FLAGS "/WX:NO ${CMAKE_EXE_LINKER_FLAGS}" ) # don't treat warnings as compile errors--gtest doesn't build
##    #set(CMAKE_CXX_FLAGS_DEBUG "/analyze ${CMAKE_CXX_FLAGS_DEBUG}")
else()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 -Wall -Werror")
    set(CMAKE_CXX_FLAGS_RELEASE "-O3")
endif()

add_subdirectory(phtree)
add_subdirectory(examples)

option(PACKAGE_TESTS "Build the tests" ON)
if(PACKAGE_TESTS)
    if (FALSE)
    
    add_compile_definitions(GTEST_HAS_ABSL=0)
    add_compile_definitions(GTEST_OS_WINDOWS_MOBILE=0)
    if(MSVC)
        add_compile_definitions(GTEST_OS_WINDOWS_MINGW=0)
    endif()
    add_compile_definitions(GTEST_OS_LINUX_ANDROID=0)
    if(LINUX) 
        add_compile_definitions(GTEST_OS_LINUX=1)
    else()
        add_compile_definitions(GTEST_OS_LINUX=0)
    endif()
    add_compile_definitions(
        GTEST_OS_WINDOWS_MOBILE=0 
        GTEST_OS_WINDOWS_PHONE=0
        GTEST_OS_WINDOWS_RT=0
        GTEST_OS_ESP8266=0
        GTEST_OS_XTENSA=0)
        
    endif()


    enable_testing()
    include(GoogleTest)
    add_subdirectory(test)
endif()
