cmake_minimum_required(VERSION 3.20)    # Minimum version to configure
project(CFD LANGUAGES CXX)   # Declare project name and language

set(CMAKE_CXX_STANDARD 20)      # Request C++ 20
set(CMAKE_CXX_STANDARD_REQUIRED ON)     # Enforce C++ 20 instead of fallback
set(CMAKE_CXX_EXTENSIONS OFF)     # Disables compiler-specific extensions 

# Create a library target named CFD
add_library(CFD
    src/mesh.cpp
    src/diffusion.cpp
    src/vtk_writer.cpp
    src/solver.cpp
)

# Add include/ to the header search path for library target CFD
# PUBLIC: The library itself can use these headers. Anything links to CFD also inherits this path
target_include_directories(CFD PUBLIC include)

# Create an executable named diffusion3d for testing
add_executable(diffusion3d apps/diffusion3d_main.cpp)
target_link_libraries(diffusion3d PRIVATE CFD)

# Create an executable named diffusion3d for simulation
add_executable(lid_driven_cavity3d apps/lid_driven_cavity3d.cpp)
target_link_libraries(lid_driven_cavity3d PRIVATE CFD)