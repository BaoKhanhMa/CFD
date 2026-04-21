#pragma once
#include "mesh.hpp"
#include "flow_state.hpp"

struct SolverSettings {
    double dt{1e-3};
    double rho{1.0};
    double nu{1e-2};
    int pressure_iters{200};
};

void apply_lid_driven_cavity_bc(const Mesh &mesh, FlowState &s);
void advance_one_step(const Mesh &mesh, const SolverSettings &cfg, FlowState &s);
double compute_div_l2(const Mesh &mesh, const FlowState &s);