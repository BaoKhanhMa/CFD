#include "vtk_writer.hpp"

#include <fstream>
#include <stdexcept>

void write_vtk(const std::string& filename, const Mesh& mesh, const FlowState& s)
{
    std::ofstream out(filename);
    if (!out) {
        throw std::runtime_error("Failed to open VTK file: " + filename);
    }

    out << "# vtk DataFile Version 3.0\n";
    out << "CFD lid driven cavity 3D\n";
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

    out << "VECTORS velocity double\n";
    for (std::size_t k = 0; k < mesh.nz(); ++k) {
        for (std::size_t j = 0; j < mesh.ny(); ++j) {
            for (std::size_t i = 0; i < mesh.nx(); ++i) {
                out << s.u(i,j,k) << ' '
                    << s.v(i,j,k) << ' '
                    << s.w(i,j,k) << '\n';
            }
        }
    }

    out << "SCALARS pressure double 1\n";
    out << "LOOKUP_TABLE default\n";
    for (std::size_t k = 0; k < mesh.nz(); ++k) {
        for (std::size_t j = 0; j < mesh.ny(); ++j) {
            for (std::size_t i = 0; i < mesh.nx(); ++i) {
                out << s.p(i,j,k) << '\n';
            }
        }
    }
}