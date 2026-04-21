#include <iostream>
#include <iomanip>
#include <sstream>

#include "mesh.hpp"
#include "flow_state.hpp"
#include "solver.hpp"
#include "vtk_writer.hpp"

int main()
{
    Mesh mesh(16, 16, 16, 1.0, 1.0, 1.0);

    FlowState s(mesh);

    SolverSettings cfg;
    cfg.dt = 1e-3;
    cfg.rho = 1.0;
    cfg.nu = 1e-2;
    cfg.pressure_iters = 30;

    apply_lid_driven_cavity_bc(mesh, s);

    write_csv_frame("cavity_0000.csv", mesh, s);

    for (int step = 1; step <= 2000; ++step) {
        advance_one_step(mesh, cfg, s);

        if (step % 20 == 0) {
            std::cout << "step=" << step
                      << " divL2=" << compute_div_l2(mesh, s)
                      << '\n';
        }

        if (step % 100 == 0) {
            std::ostringstream name;
            name << "cavity_" << std::setw(4) << std::setfill('0') << step << ".csv";
            write_csv_frame(name.str(), mesh, s);
        }
    }

    return 0;
}