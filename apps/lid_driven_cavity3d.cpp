#include <iostream>
#include "mesh.hpp"
#include "flow_state.hpp"
#include "solver.hpp"

int main(){
    Mesh mesh(24, 24, 24, 1.0, 1.0, 1.0);
    FlowState s(mesh);
    SolverSettings cfg;

    cfg.dt = 1e-3;
    cfg.rho = 1.0;
    cfg.nu = 1e-2;
    cfg.pressure_iters = 100;

    apply_lid_driven_cavity_bc(mesh, s);

    for(int step = 0; step < 5000; step++){
        advance_one_step(mesh, cfg, s);

        if(step % 20 == 0){
            std::cout << "step=" << step
                      << " divL2=" << compute_div_l2(mesh, s) << '\n';
        }
    }

    return 0;
}