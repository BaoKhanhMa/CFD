#include "mesh.hpp"

#include <vector>
#include <iostream>
#include <cmath>
#include <stdexcept>

namespace {
std::vector<double> solve_diffusion(const Mesh& mesh)
{
    const std::size_t n = mesh.cell_count();
    std::vector<double> phi(n, 0.0);
    std::vector<double> phi_new(n, 0.0);

    const double dx = mesh.dx();
    const double dy = mesh.dy();
    const double dz = mesh.dz();

    const double ax = 1.0 / (dx * dx);
    const double ay = 1.0 / (dy * dy);
    const double az = 1.0 / (dz * dz);

    const std::size_t nx = mesh.nx();
    const std::size_t ny = mesh.ny();
    const std::size_t nz = mesh.nz();

    constexpr std::size_t max_iter = 20000;
    constexpr double tol = 1e-10;

    for (std::size_t iter = 0; iter < max_iter; ++iter) {
        double max_diff = 0.0;

        for (std::size_t k = 0; k < nz; ++k) {
            for (std::size_t j = 0; j < ny; ++j) {
                for (std::size_t i = 0; i < nx; ++i) {
                    const std::size_t c = mesh.index(i, j, k);

                    if (mesh.on_xmin(i)) {
                        phi_new[c] = 1.0; // Dirichlet
                        continue;
                    }

                    if (mesh.on_xmax(i)) {
                        phi_new[c] = 0.0; // Dirichlet
                        continue;
                    }

                    const std::size_t cw = mesh.index(i - 1, j, k);
                    const std::size_t ce = mesh.index(i + 1, j, k);

                    const std::size_t cs = mesh.index(i, (j == 0 ? j : j - 1), k);
                    const std::size_t cn = mesh.index(i, (j + 1 == ny ? j : j + 1), k);

                    const std::size_t cb = mesh.index(i, j, (k == 0 ? k : k - 1));
                    const std::size_t ct = mesh.index(i, j, (k + 1 == nz ? k : k + 1));

                    const double phi_w = phi[cw];
                    const double phi_e = phi[ce];

                    // zero normal gradient on ymin/ymax/zmin/zmax
                    const double phi_s = mesh.on_ymin(j) ? phi[c] : phi[cs];
                    const double phi_n = mesh.on_ymax(j) ? phi[c] : phi[cn];
                    const double phi_b = mesh.on_zmin(k) ? phi[c] : phi[cb];
                    const double phi_t = mesh.on_zmax(k) ? phi[c] : phi[ct];

                    const double denom = 2.0 * (ax + ay + az);
                    phi_new[c] = (ax * (phi_w + phi_e)
                                + ay * (phi_s + phi_n)
                                + az * (phi_b + phi_t)) / denom;

                    max_diff = std::max(max_diff, std::abs(phi_new[c] - phi[c]));
                }
            }
        }

        phi.swap(phi_new);

        if (iter % 500 == 0) {
            std::cout << "iter=" << iter << " max_diff=" << max_diff << '\n';
        }

        if (max_diff < tol) {
            std::cout << "Converged at iter=" << iter
                      << " max_diff=" << max_diff << '\n';
            return phi;
        }
    }

    std::cout << "Warning: solver did not fully converge\n";
    return phi;
}
} // namespace

std::vector<double> run_diffusion_case(const Mesh& mesh)
{
    return solve_diffusion(mesh);
}