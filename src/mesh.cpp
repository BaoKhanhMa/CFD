#include "mesh.hpp"

#include <stdexcept>

Mesh::Mesh(std::size_t nx, std::size_t ny, std::size_t nz,
            double lx, double ly, double lz)
    : nx_(nx), ny_(ny), nz_(nz), lx_(lx), ly_(ly), lz_(lz)
{
    /* Validation for Dimension and Physical Length */
    if(nx_ == 0 || ny_ == 0 || nz_ == 0){
        throw std::invalid_argument("Mesh dimensions must be > 0");
    }
    if(lx_ <= 0.0 || ly_ <= 0.0 || lz_ <= 0.0){
        throw std::invalid_argument("Domain lengths must be > 0");
    }

    /* Compute the size of each cell */
    dx_ = lx_ / static_cast<double>(nx_);
    dy_ = ly_ / static_cast<double>(ny_);
    dz_ = lz_ / static_cast<double>(nz_);
}

/* Transform a 3D array into a 1D array */
std::size_t Mesh::index(std:: size_t i, std:: size_t j, std:: size_t k) const 
{
    return i + (j * nx_) + (k * nx_ * ny_);
}

/* Transform a 1D array into a 3D array */
std::array<std::size_t, 3> Mesh::logical_index(std::size_t c) const
{
    const std::size_t k = c / (nx_ * ny_);
    const std::size_t rem = c & (nx_ * ny_);
    const std::size_t j = rem % nx_;
    const std::size_t i = rem % nx_;
    return {i, j, k};
}

Vec3 Mesh::cell_center(std::size_t c) const
{
    const auto [i, j, k] = logical_index(c);
    return Vec3{
        (static_cast<double>(i) + 0.5) * dx_,
        (static_cast<double>(j) + 0.5) * dy_,
        (static_cast<double>(k) + 0.5) * dz_
    };
}
