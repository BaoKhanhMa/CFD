#include "vtk_writer.hpp"

#include <fstream>
#include <stdexcept>

void write_csv_frame(const std::string& filename, const Mesh& mesh, const FlowState& s)
{
    std::ofstream out(filename);
    if (!out) {
        throw std::runtime_error("Failed to open CSV file: " + filename);
    }

    out << "x,y,z,u,v,w,p\n";

    for (std::size_t k = 0; k < mesh.nz(); ++k) {
        for (std::size_t j = 0; j < mesh.ny(); ++j) {
            for (std::size_t i = 0; i < mesh.nx(); ++i) {
                const double x = (static_cast<double>(i) + 0.5) * mesh.dx();
                const double y = (static_cast<double>(j) + 0.5) * mesh.dy();
                const double z = (static_cast<double>(k) + 0.5) * mesh.dz();

                out << x << ','
                    << y << ','
                    << z << ','
                    << s.u(i,j,k) << ','
                    << s.v(i,j,k) << ','
                    << s.w(i,j,k) << ','
                    << s.p(i,j,k) << '\n';
            }
        }
    }
}