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

unit_test(test_cfilt tests/test_cfilt.c)
unit_test(test_util  tests/test_util.c)
unit_test(test_sigma tests/test_sigma.c)

binary(example_gh       examples/gh/example_gh.c)
binary(example_kalman1d examples/kalman1d/example_kalman1d.c)
binary(example_kalman   examples/kalman/example_kalman.c)
binary(example_cfilt    examples/cfilt/example_cfilt.c)
binary(example_sigma    examples/sigma/example_sigma.c)
