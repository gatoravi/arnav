cmake_minimum_required(VERSION 3.2)

project (bcall)
set (CMAKE_CXX_STANDARD 11)
include_directories ("${PROJECT_SOURCE_DIR}/src/rmath/Rmath")
add_executable (bcall "src/bcall.cc")
add_subdirectory ("${PROJECT_SOURCE_DIR}/src/")
target_link_libraries (bcall Rmath z)
