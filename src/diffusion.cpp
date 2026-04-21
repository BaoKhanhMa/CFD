#include "mesh.hpp"
#include "diffusion.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

namespace {

void apply_dirichlet_boundaries(const Mesh& mesh, std::vector<double>& phi)
{
    const std::size_t nx = mesh.nx();
    const std::size_t ny = mesh.ny();
    const std::size_t nz = mesh.nz();

    for (std::size_t k = 0; k < nz; ++k) {
        for (std::size_t j = 0; j < ny; ++j) {
            phi[mesh.index(0, j, k)] = 1.0;        // xmin
            phi[mesh.index(nx - 1, j, k)] = 0.0;   // xmax
        }
    }
}

std::vector<double> solve_diffusion(const Mesh& mesh)
{
    const std::size_t n = mesh.cell_count();
    std::vector<double> phi(n, 0.0);
    std::vector<double> phi_new(n, 0.0);

    apply_dirichlet_boundaries(mesh, phi);
    apply_dirichlet_boundaries(mesh, phi_new);

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

        apply_dirichlet_boundaries(mesh, phi_new);

        for (std::size_t k = 0; k < nz; ++k) {
            for (std::size_t j = 0; j < ny; ++j) {
                for (std::size_t i = 1; i + 1 < nx; ++i) {
                    const std::size_t c  = mesh.index(i, j, k);
                    const std::size_t cw = mesh.index(i - 1, j, k);
                    const std::size_t ce = mesh.index(i + 1, j, k);

                    const double phi_w = phi[cw];
                    const double phi_e = phi[ce];

                    const double phi_s = (j == 0)      ? phi[c] : phi[mesh.index(i, j - 1, k)];
                    const double phi_n = (j + 1 == ny) ? phi[c] : phi[mesh.index(i, j + 1, k)];
                    const double phi_b = (k == 0)      ? phi[c] : phi[mesh.index(i, j, k - 1)];
                    const double phi_t = (k + 1 == nz) ? phi[c] : phi[mesh.index(i, j, k + 1)];

                    const double denom = 2.0 * (ax + ay + az);

                    phi_new[c] = (
                        ax * (phi_w + phi_e) +
                        ay * (phi_s + phi_n) +
                        az * (phi_b + phi_t)
                    ) / denom;

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