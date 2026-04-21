#pragma once
#include "field3d.hpp"

struct FlowState {
    Field3D u;
    Field3D v;
    Field3D w;
    Field3D p;

    FlowState(const Mesh& mesh) : u(mesh, 0.0), v(mesh, 0.0), w(mesh, 0.0), p(mesh, 0.0) {};
};