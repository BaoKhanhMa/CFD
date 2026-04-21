#pragma once
#include <vector>
#include <cstddef>
#include "mesh.hpp"

class Field3D{
    public:
        Field3D() = default;
        Field3D(const Mesh& mesh, double value = 0.0) : mesh_(&mesh), values_(mesh.cell_count(), value) {}

        double& operator()(std::size_t i, std::size_t j, std::size_t k){
            return values_[mesh_->index(i, j, k)];
        }

        double operator()(std::size_t i, std::size_t j, std::size_t k) const {
            return values_[mesh_->index(i, j, k)];
        }
    
    std::vector<double>& data() { return values_; }
    const std::vector<double>& data() const { return values_; }
    
    private:
        const Mesh * mesh_{nullptr};
        std::vector<double> values_;
};