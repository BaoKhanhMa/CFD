# Physical Field

## Purpose

`Field3D` is the project’s basic scalar-field container. It stores one double-precision value for every mesh cell while letting the rest of the code read and write values using 3D indices.

Typical examples include:

- one velocity component, such as `u`
- pressure `p`
- any other cell-centered scalar quantity added later

## Design

A `Field3D` object contains two things:

- a pointer to the mesh that defines the grid layout,
- a flat `std::vector<double>` that stores the actual values.

Construction looks like this:

```cpp
Field3D(const Mesh& mesh, double value = 0.0)
```

When created, the field allocates `mesh.cell_count()` entries and initializes all of them to the requested value.

## Access pattern

The main interface is the overloaded `operator()`:

```cpp
double& operator()(std::size_t i, std::size_t j, std::size_t k);
double  operator()(std::size_t i, std::size_t j, std::size_t k) const;
```

This means a field behaves like a 3D array in user code:

```cpp
s.u(i,j,k) = 1.0;
double p = s.p(i,j,k);
```

Internally, these calls are translated into a single 1D vector access through `mesh.index(i,j,k)`.

## Internal storage

If the logical index is `(i,j,k)`, the field value lives at:

$$
q_{ijk} \longrightarrow \text{values}[i + j n_x + k n_x n_y].
$$

This gives a compact memory layout while preserving a convenient 3D programming model.

## Raw data access

`Field3D` also exposes the underlying vector through:

```cpp
std::vector<double>& data();
const std::vector<double>& data() const;
```

This is useful when:

- copying a full field,
- swapping buffers,
- writing output,
- applying lower-level numerical kernels.

## Role in the codebase

The class is intentionally minimal. It does not try to be a full tensor library or enforce physical meaning. Instead, it focuses on one job:

> Provide mesh-aware storage for a cell-centered scalar field.

That simplicity makes it easy to reuse. `FlowState` simply builds several `Field3D` instances and groups them into one fluid state object.

## Practical note

Because `Field3D` stores a pointer to an external `Mesh`, the mesh must outlive the field. In normal usage within this repository, that condition is satisfied because the mesh is created first and passed into all dependent objects.
