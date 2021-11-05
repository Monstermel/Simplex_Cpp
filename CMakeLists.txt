cmake_minimum_required(VERSION 3.22.0)
project(simplex_cpp VERSION 0.1.0)
set(CMAKE_CXX_STANDARD 17)
set (CMAKE_CXX_FLAGS "-std=c++17 -Wall -Wextra -Werror -Wshadow -pedantic")

set(Eigen_INCLUDE_DIRS "C:/Program Files/eigen-3.4.0" )
find_package(Python3 REQUIRED COMPONENTS Development)

include_directories("include/")
include_directories("${Eigen_INCLUDE_DIRS}/")
include_directories("${Python3_INCLUDE_DIRS}/")
link_directories("${Python3_LIBRARY_DIRS}/")

add_executable(solver 
    src/columnset.cpp
    src/constraint.cpp
    src/main.cpp
    src/objectivefunction.cpp
    src/simplex.cpp
    src/variable.cpp)

target_link_libraries (solver ${Python3_LIBRARIES})


