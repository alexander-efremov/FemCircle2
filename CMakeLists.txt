cmake_minimum_required(VERSION 2.8)
project(fem_circle)

set(SOURCE_FILES main.cpp)
add_executable(fem_circle ${SOURCE_FILES})

include_directories(fem)

add_subdirectory(fem)
add_subdirectory(fem_tests)
