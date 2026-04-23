# FlowState

## Purpose

`FlowState` groups together the primary variables used by the incompressible flow solver.

It contains four cell-centered fields:

- `u` : velocity in the \(x\)-direction
- `v` : velocity in the \(y\)-direction
- `w` : velocity in the \(z\)-direction
- `p` : pressure

## Structure

The type is defined as a simple struct:

```cpp
struct FlowState {
    Field3D u;
    Field3D v;
    Field3D w;
    Field3D p;
    FlowState(const Mesh& mesh);
};
```

Construction allocates all four fields on the same mesh and initializes them to zero.

## Physical meaning

Together, the velocity components form the vector field

$$
\mathbf{u} = (u,v,w).
$$

The pressure field \(p\) is used during the projection step to reduce the velocity divergence and enforce the incompressibility constraint.

In other words:

- `u`, `v`, `w` carry the momentum state,
- `p` acts as the correction field that pushes the velocity toward

$$
\nabla \cdot \mathbf{u} = 0.
$$

## Why group the variables

Bundling the flow variables into one object has a few advantages:

- solver functions can accept a single fluid-state argument,
- related fields stay synchronized on the same mesh,
- the code reads more like the physics it is implementing.

For example, the solver can naturally update a state in place:

```cpp
advance_one_step(mesh, cfg, s);
```

rather than passing four separate arrays around.

## Role in the projection method

Inside the current solver, `FlowState` is used twice:

- once for the current solution state `s`,
- once for the intermediate velocity state `star` during time advancement.

That mirrors the mathematics of a projection method:

1. compute an intermediate velocity without pressure correction,
2. solve for pressure,
3. correct the velocity field.

`FlowState` therefore acts as the numerical state vector for the CFD algorithm.
