cmake_minimum_required(VERSION 3.5)

include_directories(/home/stefan/ViennaLS/include/viennals)
include_directories(/home/stefan/ViennaHRLE/include/viennahrle)
find_package(ViennaLS REQUIRED)

add_executable("task2" task2.cpp)
