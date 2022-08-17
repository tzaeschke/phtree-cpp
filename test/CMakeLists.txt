cmake_minimum_required(VERSION 3.14)
#project(Tests)

# Avoids LNK2038 Error with MSVC
set(gtest_force_shared_crt on) 

include(FetchContent)

FetchContent_Declare(
  googletest
  GIT_REPOSITORY https://github.com/google/googletest.git
  GIT_TAG        release-1.12.1
)
FetchContent_MakeAvailable(googletest)
add_library(GTest::GTest INTERFACE IMPORTED)
target_link_libraries(GTest::GTest INTERFACE gtest_main)

#target_include_directories(${PROJECT_NAME}  PRIVATE ../phtree)


#############set(INCLUDE_DIR ${CMAKE_SOURCE_DIR}/phtree)
##############set(SOURCE_FILES phtree_test.cc ${INCLUDE_DIR}/phtree.h ${INCLUDE_DIR}/phtree_multimap.h)
############add_executable(phtree_test ${SOURCE_FILES})
###########target_include_directories(phtree_test PRIVATE ../phtree)
#target_include_directories(Tests PRIVATE ${INCLUDE_DIR})
#add_executable(phtree_test phtree_test.cc)

#target_link_libraries(phtree_test
# PRIVATE
#  GTest::GTest
#  phtree)
#############target_link_libraries(phtree_test PRIVATE GTest::GTest)


#add_test(phtree_gtests phtree_test)





# add_subdirectory("${PROJECT_SOURCE_DIR}/extern/googletest" "extern/googletest")

# The next line is optional, but keeps your CACHE cleaner:
mark_as_advanced(
    BUILD_GMOCK BUILD_GTEST BUILD_SHARED_LIBS
    gmock_build_tests gtest_build_samples gtest_build_tests
    gtest_disable_pthreads gtest_force_shared_crt gtest_hide_internal_symbols
)

# If you are interested in keeping IDEs that support folders clean, I would also add these lines:
set_target_properties(gtest PROPERTIES FOLDER extern)
set_target_properties(gtest_main PROPERTIES FOLDER extern)

macro(package_add_test TESTNAME)
    # create an executable in which the tests will be stored
    add_executable(${TESTNAME} ${ARGN})
    # link the Google test infrastructure, mocking library, and a default main function to
    # the test executable.  Remove g_test_main if writing your own main function.
    target_link_libraries(${TESTNAME} gtest gmock gtest_main)
    #target_include_directories(${TESTNAME}  PRIVATE ../phtree)
    target_include_directories(${TESTNAME} PRIVATE ..)
    # gtest_discover_tests replaces gtest_add_tests,
    # see https://cmake.org/cmake/help/v3.10/module/GoogleTest.html for more options to pass to it
    gtest_discover_tests(${TESTNAME}
        # set a working directory so your project root so that you can find test data via paths relative to the project root
        WORKING_DIRECTORY ${PROJECT_DIR}
        PROPERTIES VS_DEBUGGER_WORKING_DIRECTORY "${PROJECT_DIR}"
    )
    set_target_properties(${TESTNAME} PROPERTIES FOLDER test)
endmacro()

# package_add_test(phtree_all_test phtree_test.cc phtree_d_test.cc phtree_f_test.cc)
package_add_test(phtree_test phtree_test.cc)
package_add_test(phtree_test_const_values phtree_test_const_values.cc)
package_add_test(phtree_test_issues phtree_test_issues.cc)
package_add_test(phtree_test_ptr_values phtree_test_ptr_values.cc)
package_add_test(phtree_test_unique_ptr_values phtree_test_unique_ptr_values.cc)

package_add_test(phtree_f_test phtree_f_test.cc)

package_add_test(phtree_d_test phtree_d_test.cc)
package_add_test(phtree_d_test_copy_move phtree_d_test_copy_move.cc)
package_add_test(phtree_d_test_custom_key phtree_d_test_custom_key.cc)
package_add_test(phtree_d_test_filter phtree_d_test_filter.cc)
package_add_test(phtree_d_test_preprocessor phtree_d_test_preprocessor.cc)

package_add_test(phtree_box_f_test phtree_box_f_test.cc)

package_add_test(phtree_box_d_test phtree_box_d_test.cc)
package_add_test(phtree_box_d_test_filter phtree_box_d_test_filter.cc)
package_add_test(phtree_box_d_test_query_types phtree_box_d_test_query_types.cc)

package_add_test(phtree_multimap_d_test phtree_multimap_d_test.cc)
package_add_test(phtree_multimap_d_test_copy_move phtree_multimap_d_test_copy_move.cc)
package_add_test(phtree_multimap_d_test_filter phtree_multimap_d_test_filter.cc)
package_add_test(phtree_multimap_d_test_unique_ptr_values phtree_multimap_d_test_unique_ptr_values.cc)

package_add_test(phtree_multimap_box_d_test phtree_multimap_box_d_test.cc)
