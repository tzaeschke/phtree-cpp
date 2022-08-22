cmake_minimum_required(VERSION 3.14)
project(phtree-tests LANGUAGES CXX)

# Avoids LNK2038 Error with MSVC
set(gtest_force_shared_crt on)

include(FetchContent)
include(common/scripts.cmake)

FetchContent_Declare(
        googletest
        GIT_REPOSITORY https://github.com/google/googletest.git
        GIT_TAG release-1.12.1
)
FetchContent_MakeAvailable(googletest)
add_library(GTest::GTest INTERFACE IMPORTED)
target_link_libraries(GTest::GTest INTERFACE gtest_main)

# The next line is optional, but keeps your CACHE cleaner:
mark_as_advanced(
        BUILD_GMOCK BUILD_GTEST BUILD_SHARED_LIBS
        gmock_build_tests gtest_build_samples gtest_build_tests
        gtest_disable_pthreads gtest_force_shared_crt gtest_hide_internal_symbols
)

# If you are interested in keeping IDEs that support folders clean, I would also add these lines:
set_target_properties(gtest PROPERTIES FOLDER extern)
set_target_properties(gtest_main PROPERTIES FOLDER extern)


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

add_subdirectory(common)

