cmake_minimum_required(VERSION 2.8)

project(fem)

add_definitions(-std=c++0x)


set(HEADER_FILES ${CMAKE_CURRENT_SOURCE_DIR})

find_package( OpenMP REQUIRED)
if(OPENMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

set(BOOST_ROOT "C:\\boost_1_62_0")

find_package(Boost 1.62.0)
if(Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS})
endif()

set(SOURCE_FILES graphs.cpp solver1.cpp solver2.cpp consts.cpp)

add_library(fem STATIC ${SOURCE_FILES} ${HEADER_FILES})


if(OPENMP_FOUND)
    target_link_libraries(fem ${OpenMP_CXX_FLAGS})
endif()

if(Boost_FOUND)
    target_link_libraries(fem ${Boost_LIBRARIES})
endif()


