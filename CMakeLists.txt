cmake_minimum_required(VERSION 3.14.5)

project(standalone_operators LANGUAGES CXX)

include(FetchContent)
FetchContent_Declare(
    range-v3
    GIT_REPOSITORY https://github.com/ericniebler/range-v3.git
    GIT_TAG         0.5.0
    )
FetchContent_MakeAvailable(range-v3)

add_executable(standalone main.cc)
target_compile_features(standalone PRIVATE cxx_std_17)
target_link_libraries(standalone PRIVATE range-v3)

add_executable(test_center test_center.cc)
target_compile_features(test_center PRIVATE cxx_std_17)
target_link_libraries(test_center PRIVATE range-v3)

add_executable(diff_explicit diff_explicit.cc)
target_compile_features(diff_explicit PRIVATE cxx_std_17)
target_link_libraries(diff_explicit PRIVATE range-v3)

add_executable(test_ranges test_ranges.cc)
target_compile_features(test_ranges PRIVATE cxx_std_17)
target_link_libraries(test_ranges PRIVATE range-v3)

set(bla ${CMAKE_CURRENT_BINARY_DIR})
