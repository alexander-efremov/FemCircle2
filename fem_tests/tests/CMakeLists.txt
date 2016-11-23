include_directories(${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})

add_definitions(-std=c++0x)

add_executable(runBasicFemTests
        tests.cpp)

FIND_PACKAGE( OpenMP REQUIRED)
if(OPENMP_FOUND)
    message("OPENMP FOUND")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

set(BOOST_ROOT "C:\\boost_1_62_0")

find_package(Boost 1.62.0 COMPONENTS system filesystem REQUIRED)
if(Boost_FOUND)
    message(STATUS "Boost_INCLUDE_DIRS: ${Boost_INCLUDE_DIRS}")
    message(STATUS "Boost_LIBRARIES: ${Boost_LIBRARIES}")
    message(STATUS "Boost_VERSION: ${Boost_VERSION}")
    include_directories(${Boost_INCLUDE_DIRS})
endif()

target_link_libraries(runBasicFemTests gtest gtest_main)
if(Boost_FOUND)
    target_link_libraries(runBasicFemTests ${Boost_LIBRARIES})
endif()
if(OPENMP_FOUND)
    target_link_libraries(runBasicFemTests fem ${OpenMP_CXX_FLAGS})
endif()