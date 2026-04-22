

#pragma once

#include <array>
#include <vector>
#include <cstddef>

/* A simple vector struct with value-initialization as 0.0 by default */
struct Vec3{
    double x{};
    double y{};
    double z{};
};

class Mesh {
    public:
        /* Prototype for constructor */
        Mesh(std::size_t nx, std::size_t ny, std::size_t nz,
                double lx, double ly, double lz);

        std::size_t nx() const { return nx_; }
        std::size_t ny() const { return ny_; }
        std::size_t nz() const { return nz_; }

        std::size_t cell_count() const { return nx_ * ny_ * nz_;}
        
        double dx() const { return dx_; }
        double dy() const { return dy_; }
        double dz() const { return dz_; }

        double cell_volume() const { return dx_ * dy_ * dz_; }

        /* Prototypes */
        std::size_t index(std::size_t i, std::size_t j, std::size_t k) const;
        std::array<std::size_t, 3> logical_index(std::size_t c) const;

        Vec3 cell_center(std::size_t c) const;

        bool on_xmin(std::size_t i) const { return i == 0; }
        bool on_xmax(std::size_t i) const { return i + 1 == nx_; }
        bool on_ymin(std::size_t j) const { return j == 0; }
        bool on_ymax(std::size_t j) const { return j + 1 == ny_; }
        bool on_zmin(std::size_t k) const { return k == 0; }
        bool on_zmax(std::size_t k) const { return k + 1 == nz_; }

    private:
        /* How many cells in each dimension? */
        std::size_t nx_{};
        std::size_t ny_{};
        std::size_t nz_{};

        /* What is the physical size in each dimension? */
        double lx_{};
        double ly_{};
        double lz_{};

        /* What is the size of each cell in each dimension? */
        double dx_{};
        double dy_{};
        double dz_{};
};