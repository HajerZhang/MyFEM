cmake_minimum_required(VERSION 3.12)

project(
    MyFEM
    VERSION 0.1
    LANGUAGES CXX
)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(srcs)
add_subdirectory("src")
include_directories("include")

set(CMAKE_PREFIX_PATH "/opt/homebrew/share/eigen3/cmake")
find_package(Eigen3)
if(Eigen3_FOUND)
    message(STATUS "Eigen3 found: version ${Eigen3_VERSION}") 
    include_directories(${EIGEN3_INCLUDE_DIR})
else()
    message(FATAL_ERROR "Eigen3 not found")
endif()

add_executable(MyFEM ${srcs})
