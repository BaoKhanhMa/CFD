#include "mesh.hpp"
#include "diffusion.hpp"

#include <iostream>
#include <vector>
#include <fstream>
#include <stdexcept>

static void write_legacy_vtk(const Mesh& mesh, const std::vector<double>& phi)
{
    std::ofstream out("diffusion3d.vtk");
    if (!out) {
        throw std::runtime_error("Failed to open diffusion3d.vtk");
    }

    out << "# vtk DataFile Version 3.0\n";
    out << "3D diffusion result\n";
    out << "ASCII\n";
    out << "DATASET STRUCTURED_POINTS\n";
    out << "DIMENSIONS " << mesh.nx() + 1 << ' '
                        << mesh.ny() + 1 << ' '
                        << mesh.nz() + 1 << '\n';
    out << "ORIGIN 0 0 0\n";
    out << "SPACING " << mesh.dx() << ' '
                      << mesh.dy() << ' '
                      << mesh.dz() << '\n';
    out << "CELL_DATA " << mesh.cell_count() << '\n';
    out << "SCALARS phi double\n";
    out << "LOOKUP_TABLE default\n";

    for (double v : phi) {
        out << v << '\n';
    }
}

int main()
{
    try {
        Mesh mesh(32, 32, 32, 1.0, 1.0, 1.0);

        const auto phi = run_diffusion_case(mesh);

        write_legacy_vtk(mesh, phi);

        std::cout << "Wrote diffusion3d.vtk\n";
        return 0;
    }
    catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << '\n';
        return 1;
    }
}