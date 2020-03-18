cmake_minimum_required (VERSION 3.11)
cmake_policy(VERSION 3.11)

find_package(Catch2 REQUIRED)

set(GOMA_TEST_SOURCES
    gds/gds_vector.cpp
)
add_executable(goma_unit_tests tests.cpp ${GOMA_TEST_SOURCES})
target_link_libraries(goma_unit_tests Catch2::Catch2 gds)

include(CTest)
include(Catch)
catch_discover_tests(goma_unit_tests)