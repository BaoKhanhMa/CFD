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

    private:
        std::size_t nx_{};
        std::size_t ny_{};
        std::size_t nz_{};

        double lx_{};
        double ly_{};
        double lz_{};

        double dx_{};
        double dy_{};
        double dz_{};
};