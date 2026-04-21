#pragma once

#include <string>

#include "mesh.hpp"
#include "flow_state.hpp"

void write_vtk(const std::string& filename, const Mesh& mesh, const FlowState& s);