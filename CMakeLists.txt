cmake_minimum_required(VERSION 3.5)

include_directories(/home/stefan/ViennaLS/include/viennals)
include_directories(/home/stefan/ViennaHRLE/include/viennahrle)
find_package(ViennaLS REQUIRED)

add_executable("frosty" Frosty.cpp)
add_executable("airgap" airgap.cpp)
add_executable("deposition" deposition.cpp)
add_executable("geomadvection" geomadvection.cpp)
add_executable("periodBC" periodBC.cpp)
add_executable("test" Test.cpp)
add_executable("pattern" pattern.cpp)
