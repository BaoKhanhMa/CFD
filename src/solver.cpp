#include "solver.hpp"
#include <cmath>
#include <algorithm>

static void copy_field(const std::vector<double>& src, std::vector<double>& dst) {
    dst = src;
}

void apply_lid_driven_cavity_bc(const Mesh& mesh, FlowState& s) {
    for (std::size_t k = 0; k < mesh.nz(); ++k) {
        for (std::size_t j = 0; j < mesh.ny(); ++j) {
            for (std::size_t i = 0; i < mesh.nx(); ++i) {
                if (mesh.on_xmin(i) || mesh.on_xmax(i) ||
                    mesh.on_ymin(j) || mesh.on_ymax(j) ||
                    mesh.on_zmin(k)) {
                    s.u(i,j,k) = 0.0;
                    s.v(i,j,k) = 0.0;
                    s.w(i,j,k) = 0.0;
                }
                if (mesh.on_zmax(k)) {
                    s.u(i,j,k) = 1.0; // moving lid
                    s.v(i,j,k) = 0.0;
                    s.w(i,j,k) = 0.0;
                }
            }
        }
    }
}

static double laplacian(const Field3D& f, const Mesh& mesh,
                        std::size_t i, std::size_t j, std::size_t k) {
    const double dx2 = mesh.dx() * mesh.dx();
    const double dy2 = mesh.dy() * mesh.dy();
    const double dz2 = mesh.dz() * mesh.dz();

    const double xc = f(i,j,k);
    const double xm = (i > 0) ? f(i-1,j,k) : xc;
    const double xp = (i+1 < mesh.nx()) ? f(i+1,j,k) : xc;
    const double ym = (j > 0) ? f(i,j-1,k) : xc;
    const double yp = (j+1 < mesh.ny()) ? f(i,j+1,k) : xc;
    const double zm = (k > 0) ? f(i,j,k-1) : xc;
    const double zp = (k+1 < mesh.nz()) ? f(i,j,k+1) : xc;

    return (xp - 2.0*xc + xm)/dx2
         + (yp - 2.0*xc + ym)/dy2
         + (zp - 2.0*xc + zm)/dz2;
}

static double divergence(const FlowState& s, const Mesh& mesh,
                         std::size_t i, std::size_t j, std::size_t k) {
    const double um = (i > 0) ? s.u(i-1,j,k) : s.u(i,j,k);
    const double up = (i+1 < mesh.nx()) ? s.u(i+1,j,k) : s.u(i,j,k);

    const double vm = (j > 0) ? s.v(i,j-1,k) : s.v(i,j,k);
    const double vp = (j+1 < mesh.ny()) ? s.v(i,j+1,k) : s.v(i,j,k);

    const double wm = (k > 0) ? s.w(i,j,k-1) : s.w(i,j,k);
    const double wp = (k+1 < mesh.nz()) ? s.w(i,j,k+1) : s.w(i,j,k);

    return (up - um) / (2.0 * mesh.dx())
         + (vp - vm) / (2.0 * mesh.dy())
         + (wp - wm) / (2.0 * mesh.dz());
}

static void pressure_poisson(const Mesh& mesh, const SolverSettings& cfg,
                             const FlowState& star, Field3D& p) {
    Field3D p_new(mesh, 0.0);

    const double dx2 = mesh.dx() * mesh.dx();
    const double dy2 = mesh.dy() * mesh.dy();
    const double dz2 = mesh.dz() * mesh.dz();
    const double denom = 2.0/dx2 + 2.0/dy2 + 2.0/dz2;

    for (int it = 0; it < cfg.pressure_iters; ++it) {
        for (std::size_t k = 0; k < mesh.nz(); ++k) {
            for (std::size_t j = 0; j < mesh.ny(); ++j) {
                for (std::size_t i = 0; i < mesh.nx(); ++i) {
                    const double pc = p(i,j,k);
                    const double pmx = (i > 0) ? p(i-1,j,k) : pc;
                    const double ppx = (i+1 < mesh.nx()) ? p(i+1,j,k) : pc;
                    const double pmy = (j > 0) ? p(i,j-1,k) : pc;
                    const double ppy = (j+1 < mesh.ny()) ? p(i,j+1,k) : pc;
                    const double pmz = (k > 0) ? p(i,j,k-1) : pc;
                    const double ppz = (k+1 < mesh.nz()) ? p(i,j,k+1) : pc;

                    const double rhs = (cfg.rho / cfg.dt) * divergence(star, mesh, i, j, k);

                    p_new(i,j,k) = (
                        (pmx + ppx)/dx2 +
                        (pmy + ppy)/dy2 +
                        (pmz + ppz)/dz2 - rhs
                    ) / denom;
                }
            }
        }
        
        p.data().swap(p_new.data());
    }
}

void advance_one_step(const Mesh& mesh, const SolverSettings& cfg, FlowState& s) {
    FlowState star(mesh);

    for (std::size_t k = 0; k < mesh.nz(); ++k) {
        for (std::size_t j = 0; j < mesh.ny(); ++j) {
            for (std::size_t i = 0; i < mesh.nx(); ++i) {
                star.u(i,j,k) = s.u(i,j,k) + cfg.dt * cfg.nu * laplacian(s.u, mesh, i,j,k);
                star.v(i,j,k) = s.v(i,j,k) + cfg.dt * cfg.nu * laplacian(s.v, mesh, i,j,k);
                star.w(i,j,k) = s.w(i,j,k) + cfg.dt * cfg.nu * laplacian(s.w, mesh, i,j,k);
            }
        }
    }

    apply_lid_driven_cavity_bc(mesh, star);
    pressure_poisson(mesh, cfg, star, s.p);

    for (std::size_t k = 0; k < mesh.nz(); ++k) {
        for (std::size_t j = 0; j < mesh.ny(); ++j) {
            for (std::size_t i = 0; i < mesh.nx(); ++i) {
                const double pmx = (i > 0) ? s.p(i-1,j,k) : s.p(i,j,k);
                const double ppx = (i+1 < mesh.nx()) ? s.p(i+1,j,k) : s.p(i,j,k);
                const double pmy = (j > 0) ? s.p(i,j-1,k) : s.p(i,j,k);
                const double ppy = (j+1 < mesh.ny()) ? s.p(i,j+1,k) : s.p(i,j,k);
                const double pmz = (k > 0) ? s.p(i,j,k-1) : s.p(i,j,k);
                const double ppz = (k+1 < mesh.nz()) ? s.p(i,j,k+1) : s.p(i,j,k);

                s.u(i,j,k) = star.u(i,j,k) - cfg.dt/cfg.rho * (ppx - pmx) / (2.0 * mesh.dx());
                s.v(i,j,k) = star.v(i,j,k) - cfg.dt/cfg.rho * (ppy - pmy) / (2.0 * mesh.dy());
                s.w(i,j,k) = star.w(i,j,k) - cfg.dt/cfg.rho * (ppz - pmz) / (2.0 * mesh.dz());
            }
        }
    }

    apply_lid_driven_cavity_bc(mesh, s);
}

double compute_div_l2(const Mesh& mesh, const FlowState& s) {
    double sum = 0.0;
    for (std::size_t k = 0; k < mesh.nz(); ++k) {
        for (std::size_t j = 0; j < mesh.ny(); ++j) {
            for (std::size_t i = 0; i < mesh.nx(); ++i) {
                const double d = divergence(s, mesh, i, j, k);
                sum += d * d;
            }
        }
    }
    return std::sqrt(sum / static_cast<double>(mesh.cell_count()));
}