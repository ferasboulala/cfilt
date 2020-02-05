cmake_minimum_required(VERSION 2.8)
# Using the most recent standard
set (CMAKE_C_FLAGS "-std=gnu99 ${CMAKE_C_FLAGS}")

# Warnings as errors
set(CMAKE_C_FLAGS "-Werror -Wall ${CMAKE_C_FLAGS}")

# Release build
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS_RELEASE} ${CMAKE_C_FLAGS}")

# Outputing binaries in bin directory
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/bin)

project(cfilt C)

include_directories(. examples tests)
file(GLOB SRC_FILES cfilt/*.c)

macro(binary prog_name prog_src)
    add_executable(${prog_name} ${prog_src} ${SRC_FILES})
    target_link_libraries(${prog_name} m gsl gslcblas)
endmacro()

enable_testing()
macro(unit_test test_name test_source)
    add_executable(${test_name} ${test_source} ${SRC_FILES})
    add_test(${test_name} ${CMAKE_CURRENT_SOURCE_DIR}/bin/${test_name})
    target_link_libraries(${test_name} m gsl gslcblas)
endmacro()

unit_test(test_cfilt  tests/test_cfilt.c)
unit_test(test_util   tests/test_util.c)
unit_test(test_sigma  tests/test_sigma.c)
unit_test(test_gh     tests/test_gh.c)
unit_test(test_kalman tests/test_kalman.c)
unit_test(test_ukf    tests/test_ukf.c)

binary(discrete_white_noise examples/cfilt/discrete_white_noise.c)
binary(mahalanobis          examples/cfilt/mahalanobis.c)
binary(nees                 examples/cfilt/nees.c)
binary(van_der_merwe        examples/sigma/van_der_merwe.c)
