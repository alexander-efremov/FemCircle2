cmake_minimum_required(VERSION 2.8)

project(fem)

add_definitions(-std=c++0x)

set(HEADER_FILES
    utils.h
    common.h
    consts.h
    )

FIND_PACKAGE( OpenMP REQUIRED)
if(OPENMP_FOUND)
    message("OPENMP FOUND")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

set(SOURCE_FILES
        solver1.cpp solver2.cpp consts.cpp)

add_library(fem STATIC ${SOURCE_FILES} ${HEADER_FILES})
target_link_libraries(fem ${OpenMP_CXX_FLAGS})
