cmake_minimum_required(VERSION 3.21)
project(Elliptic_SFA)

set(CMAKE_CXX_STANDARD 14)

include_directories(.)
find_library(libconfig libconfig)

add_executable(Elliptic_SFA
        EllipSFA.cpp
        EllipSFA.hpp
        main.cpp)

target_link_libraries(Elliptic_SFA config++)