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

include_directories(. tests)
file(GLOB SRC_FILES cfilt/*.c)

macro(binary prog_name prog_src)
    add_executable(${prog_name} ${prog_src} ${SRC_FILES})
    target_link_libraries(${prog_name} m gsl gslcblas)
endmacro()

binary(test_gh tests/test_gh.c)
binary(test_kalman1d tests/test_kalman1d.c)
binary(test_kalman tests/test_kalman.c)
