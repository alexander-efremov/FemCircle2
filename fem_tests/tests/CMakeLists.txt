include_directories(${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})

add_definitions(-std=c++0x)


if (CMAKE_HOST_WIN32)
    set(BOOST_ROOT "C:\\boost_1_62_0")
else()
    set(BOOST_ROOT "/home/jane/ClionProjects/boost_1_62_0")
endif()

set(CPP_SOURCES tests.cpp graph_tests.cpp graph_tests.cpp)
add_executable(runBasicFemTests ${CPP_SOURCES})

find_package( OpenMP REQUIRED)
if(OPENMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

find_package(Boost 1.62.0)
if(Boost_FOUND)
    message(${Boost_INCLUDE_DIRS})
    include_directories(${Boost_INCLUDE_DIRS})
endif()

target_link_libraries(runBasicFemTests gtest gtest_main)

if(Boost_FOUND)
    target_link_libraries(runBasicFemTests ${Boost_LIBRARIES})
endif()
if(OPENMP_FOUND)
    target_link_libraries(runBasicFemTests fem ${OpenMP_CXX_FLAGS})
endif()


