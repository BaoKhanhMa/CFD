# Mesh manual

## Purpose

`Mesh` defines a uniform 3D Cartesian grid over a rectangular domain. It stores:

- the number of cells in each direction: `nx`, `ny`, `nz`
- the physical domain lengths: `lx`, `ly`, `lz`
- the derived cell sizes: `dx`, `dy`, `dz`

This class is the geometric foundation for the rest of the code. Every field and solver routine uses the mesh to interpret array indices as physical locations.

## Domain and spacing

The constructor takes:

```cpp
Mesh(std::size_t nx, std::size_t ny, std::size_t nz,
     double lx, double ly, double lz);
```

From these values, the mesh computes uniform spacing:

$$
\Delta x = \frac{L_x}{n_x}, \qquad
\Delta y = \frac{L_y}{n_y}, \qquad
\Delta z = \frac{L_z}{n_z}.
$$

Each cell therefore has volume:

$$
V_{cell} = (\Delta x \, \Delta y \, \Delta z).
$$

The constructor also validates that all cell counts are positive and that all physical lengths are strictly positive.

## Storage model

Although the mesh is logically 3D, the project stores field values in a 1D array. `Mesh` provides the mapping between the two views.

### 3D to 1D indexing

A cell at logical coordinates `(i,j,k)` is flattened as:

$$
\text{index}(i, j, k) = i + (j\*n_x) + (k\*n_x*n_y).
$$

Interpretation:

- moving one cell in `x` increases the index by `1`
- moving one cell in `y` increases the index by `n_x`
- moving one cell in `z` increases the index by `n_x n_y`

This is a standard row-major layout for a structured grid.

## 1D to 3D indexing

The inverse map is implemented by `logical_index(c)`, which reconstructs `(i,j,k)` from a flat cell index `c`.

Conceptually:

$$
k = \left\lfloor \frac{c}{n_x n_y} \right\rfloor
$$

$$
\text{remainder} = c \bmod (n_x n_y)
$$

$$
j = \left\lfloor \frac{\text{remainder}}{n_x} \right\rfloor
$$

$$
\qquad
i = \text{remainder} \bmod n_x.
$$

This is useful when iterating linearly over stored data but still needing geometric coordinates.

## Cell centers

`cell_center(c)` returns the physical center of cell `c` as a `Vec3`:

$$
x_c = \left(i + \tfrac{1}{2}\right)\Delta x,
$$
$$
\qquad
y_c = \left(j + \tfrac{1}{2}\right)\Delta y,
$$
$$
\qquad
z_c = \left(k + \tfrac{1}{2}\right)\Delta z.
$$

This is the natural location for cell-centered variables such as pressure and the velocity components in this codebase.

## Boundary helpers

The mesh also exposes simple boundary tests:

- `on_xmin(i)`, `on_xmax(i)`
- `on_ymin(j)`, `on_ymax(j)`
- `on_zmin(k)`, `on_zmax(k)`

These helpers keep boundary-condition code readable and make it explicit when a solver step is acting on a wall or lid.

## Why this class matters

`Mesh` is deliberately small, but it carries three essential responsibilities:

1. it defines the computational domain,
2. it gives every field a consistent indexing convention,
3. it supplies geometric factors needed by finite-difference operators.

Without this layer, every solver routine would need to re-implement indexing, spacing, and boundary detection on its own.
