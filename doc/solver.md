# Solver

## Definition

The solver advances a 3D incompressible lid-driven cavity flow on a uniform Cartesian mesh. The implementation follows a simple projection-style method:

1. apply velocity boundary conditions,
2. compute an intermediate velocity field,
3. solve a Poisson equation for pressure,
4. correct the velocity using the pressure gradient,
5. re-apply boundary conditions.

The code is intentionally compact, so the numerical method is easier to read than optimize.

## Governing equations

The target physics is the incompressible Navier-Stokes system:

$$
\nabla \cdot \mathbf{u} = 0,
$$

$$
\frac{\partial \mathbf{u}}{\partial t} + (\mathbf{u}\cdot\nabla)\mathbf{u}
= -\frac{1}{\rho}\nabla p + \nu \nabla^2 \mathbf{u}.
$$

Here:

- $\mathbf{u} = (u,v,w)$ is velocity,
- $p$ is pressure,
- $\rho$ is density,
- $\nu$ is kinematic viscosity.

The solver settings expose the main numerical and physical parameters:

- `dt` : time step $\Delta t$
- `rho` : density $\rho$
- `nu` : kinematic viscosity $\nu$
- `pressure_iters` : number of Poisson iterations per time step

## Boundary conditions

The implemented case is a **lid-driven cavity**.

Five faces are stationary no-slip walls:

$$
\mathbf{u} = (0,0,0).
$$

The top face $(z = L_z)$ is a moving lid with unit tangential velocity:

$$
u = 1, \qquad v = 0, \qquad w = 0.
$$

This boundary condition injects momentum into the cavity and generates the recirculating flow.

## Time advancement

### 1. Intermediate velocity

The code first computes an intermediate state $\mathbf{u}^*$ without the pressure-gradient term:

$$
\frac{\mathbf{u}^* - \mathbf{u}^n}{\Delta t}
= - (\mathbf{u}^n \cdot \nabla)\mathbf{u}^n + \nu \nabla^2 \mathbf{u}^n.
$$

Equivalently,

$$
\mathbf{u}^* = \mathbf{u}^n + \Delta t
\left[- (\mathbf{u}^n \cdot \nabla)\mathbf{u}^n + \nu \nabla^2 \mathbf{u}^n\right].
$$

The same update is applied component-wise to `u`, `v`, and `w`.

### 2. Pressure equation

To enforce incompressibility after correction, pressure is obtained from a Poisson equation driven by the divergence of the intermediate field:

$$
\nabla^2 p = \frac{\rho}{\Delta t} \nabla \cdot \mathbf{u}^*.
$$

The code solves this equation iteratively with repeated stencil sweeps.

### 3. Velocity correction

Once pressure is available, the velocity is projected as:

$$
\mathbf{u}^{n+1} = \mathbf{u}^* - \frac{\Delta t}{\rho} \nabla p.
$$

This step removes much of the divergence from the intermediate field and moves the state toward

$$
\nabla \cdot \mathbf{u}^{n+1} = 0.
$$

## Discrete operators

All operators are evaluated on a uniform cell-centered mesh.

### Laplacian

For a scalar field $\phi$, the code uses the standard second-order central form:

$$
\nabla^2 \phi \approx
\frac{\phi_{i+1,j,k} - 2\phi_{i,j,k} + \phi_{i-1,j,k}}{\Delta x^2}
+
\frac{\phi_{i,j+1,k} - 2\phi_{i,j,k} + \phi_{i,j-1,k}}{\Delta y^2}
+
\frac{\phi_{i,j,k+1} - 2\phi_{i,j,k} + \phi_{i,j,k-1}}{\Delta z^2}.
$$

At domain boundaries, the implementation falls back to the center value when a neighbor does not exist. In practice, that acts like a simple one-sided closure tied to the imposed boundary value.

### Advection

The convective term is evaluated component-wise. For example, the `u`-equation uses

$$
(\mathbf{u}\cdot\nabla)u = u\frac{\partial u}{\partial x} + v\frac{\partial u}{\partial y} + w\frac{\partial u}{\partial z}.
$$

The directional derivative in each coordinate is approximated with a first-order **upwind** difference chosen by the sign of the transporting velocity. In 1D notation,

$$
\frac{\partial \phi}{\partial x} \approx
\begin{cases}
\dfrac{\phi_i - \phi_{i-1}}{\Delta x}, & u \ge 0, \\
\dfrac{\phi_{i+1} - \phi_i}{\Delta x}, & u < 0.
\end{cases}
$$

This is more dissipative than a centered scheme, but it is robust and simple for an experimental solver.

### Divergence

The divergence monitor uses central differences:

$$
\nabla \cdot \mathbf{u} \approx
\frac{u_{i+1,j,k} - u_{i-1,j,k}}{2\Delta x}
+
\frac{v_{i,j+1,k} - v_{i,j-1,k}}{2\Delta y}
+
\frac{w_{i,j,k+1} - w_{i,j,k-1}}{2\Delta z}.
$$

The helper `compute_div_l2()` reports the root-mean-square divergence over all cells:

$$
\|\nabla \cdot \mathbf{u}\|_{L_2}
\approx
\sqrt{\frac{1}{N}\sum_{c=1}^{N}\left(\nabla \cdot \mathbf{u}\right)_c^2 }.
$$

This gives a simple measure of how well incompressibility is being enforced.

## Pressure Poisson iteration

The Poisson equation is solved with repeated explicit stencil updates. For each cell,

```math
p_{i,j,k}^{\mathrm{new}} =
\frac{
\frac{p_{i-1,j,k}+p_{i+1,j,k}}{\Delta x^2} +
\frac{p_{i,j-1,k}+p_{i,j+1,k}}{\Delta y^2} +
\frac{p_{i,j,k-1}+p_{i,j,k+1}}{\Delta z^2}
-
\left(\frac{\rho}{\Delta t}\nabla \cdot \mathbf{u}^*\right)_{i,j,k}
}{
\frac{2}{\Delta x^2} + \frac{2}{\Delta y^2} + \frac{2}{\Delta z^2}
}
```

This is a simple Jacobi-style relaxation. The number of iterations is fixed by `pressure_iters`.

## Algorithm summary

At each time step, the implemented routine `advance_one_step()` performs:

1. build a temporary flow state `star`,
2. update `star.u`, `star.v`, `star.w` with advection and diffusion,
3. apply cavity boundary conditions to `star`,
4. solve the pressure Poisson equation into `s.p`,
5. correct `s.u`, `s.v`, `s.w` using the pressure gradient,
6. re-apply cavity boundary conditions to the corrected state.

In compact notation:

$$
\mathbf{u}^* \leftarrow \text{advect-diffuse}(\mathbf{u}^n),
$$

$$
p^{n+1} \leftarrow \text{solve Poisson}(\nabla\cdot\mathbf{u}^*),
$$

$$
\mathbf{u}^{n+1} \leftarrow \mathbf{u}^* - \frac{\Delta t}{\rho}\nabla p^{n+1}.
$$

## What this solver is and is not

This solver is a clear educational baseline. It already contains the essential CFD ideas:

- structured mesh storage,
- finite-difference operators,
- explicit advection-diffusion update,
- projection-based pressure correction,
- physically meaningful cavity boundary conditions.

At the same time, it is still a simplified implementation. For example:

- all variables are cell-centered,
- the pressure solve uses a fixed iteration count instead of a true residual-based convergence check,
- the advection scheme is first-order upwind,
- boundary handling is simple rather than high-order.

That makes it a good platform for experimentation, extension, and learning.
