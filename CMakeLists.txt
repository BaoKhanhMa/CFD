cmake_minimum_required(VERSION 3.20)   
project(CFD LANGUAGES CXX) 

set(CMAKE_CXX_STANDARD 20)     
set(CMAKE_CXX_STANDARD_REQUIRED ON)    
set(CMAKE_CXX_EXTENSIONS OFF)   

add_library(CFD
    src/mesh.cpp
    src/solver.cpp
)

target_include_directories(CFD PUBLIC include)

add_executable(lid_driven_cavity3d apps/lid_driven_cavity3d.cpp)
target_link_libraries(lid_driven_cavity3d PRIVATE CFD)