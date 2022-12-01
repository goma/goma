cmake_minimum_required (VERSION 3.11)
cmake_policy(VERSION 3.11)

find_package(Catch2 3 REQUIRED)

set(GOMA_TEST_SOURCES
    gds/gds_vector.cpp
    bc/rotate_util.cpp
)

add_executable(goma_unit_tests unit_tests_main.cpp ${GOMA_TEST_SOURCES})
target_link_libraries(goma_unit_tests Catch2::Catch2 goma_util gds ${MPI_C_LIBRARIES})
target_include_directories(goma_unit_tests PRIVATE ${MPI_C_INCLUDE_PATH})

include(CTest)
include(Catch)
catch_discover_tests(goma_unit_tests)
