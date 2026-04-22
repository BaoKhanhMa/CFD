# CFD

## Introduction

This project is a small 3D CFD code in C++ built on a structured Cartesian mesh.

The current simulation path is a transient lid-driven cavity solver advanced in time step by step. The solver stores velocity and pressure at cell centers and uses a projection-style update to enforce incompressibility.

## System Design

The code is split into a few small modules. Each module has its own manual.

- [Mesh](doc/mesh.md)  
  Defines the structured 3D grid, geometric spacing, flat indexing, cell centers, and boundary checks.

- [Field3D](doc/field3d.md)  
  Defines scalar storage on the mesh and the memory layout used by the solver.

- [FlowState](doc/flow_state.md)  
  Groups the flow unknowns `u`, `v`, `w`, and `p` into one state object.

- [Solver](doc/solver.md)  
  Defines the time-marching algorithm, lid-driven cavity boundary conditions, pressure projection, and divergence monitor.

## Build

```bash
mkdir -p build
cd build
cmake ..
cmake --build .
