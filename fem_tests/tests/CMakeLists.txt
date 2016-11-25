include_directories(${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})

add_definitions(-std=c++0x)


set(BOOST_ROOT "C:\\boost_1_62_0" )

link_directories(${GVZ_LIB_PATH})

set(CPP_SOURCES tests.cpp)
set(PY_SOURCES graphviz.py)
add_executable(runBasicFemTests ${CPP_SOURCES} ${PY_SOURCES})

find_package( OpenMP REQUIRED)
if(OPENMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

find_package(Boost 1.62.0 )
if(Boost_FOUND)
    message(${Boost_INCLUDE_DIRS})
    include_directories(${Boost_INCLUDE_DIRS})
endif()

#find_library(GRAPVIZ_GVC gvc HINTS ${GVZ_LIB_PATH})
#find_library(GRAPVIZ_CGRAPH cgraph HINTS ${GVZ_LIB_PATH})
#find_library(GRAPVIZ_CDT cdt HINTS ${GVZ_LIB_PATH})
#message(WARNING "GRAPVIZ_GVC is: ${GRAPVIZ_GVC}")
#message(WARNING "GRAPVIZ_CGRAPH is: ${GRAPVIZ_CGRAPH}")
#message(WARNING "GRAPVIZ_CDT is: ${GRAPVIZ_CDT}")

target_link_libraries(runBasicFemTests gtest gtest_main )
#target_link_libraries(runBasicFemTests gvc )

if(Boost_FOUND)
    target_link_libraries(runBasicFemTests ${Boost_LIBRARIES})
endif()
if(OPENMP_FOUND)
    target_link_libraries(runBasicFemTests fem ${OpenMP_CXX_FLAGS})
endif()


