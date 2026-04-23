# Diffusion manual

## Problem being solved

The diffusion example solves a steady scalar problem in a 3D box. The unknown scalar field $\phi$ can be interpreted as temperature, concentration, or electric potential.

The current setup uses:

- left face $(x = 0)$:

  $$\phi = 1$$

- right face $(x = L_x)$:

  $$\phi = 0$$

- other faces: zero normal flux, represented in the implementation by reusing the cell value when a neighbor outside the domain would be needed.

## Governing equation

The code solves the Laplace equation:

$$
\nabla^2 \phi = 0
$$

This is the steady diffusion equation with no source term.

## Discretization

On the uniform mesh, the discrete equation at each interior cell is:

$$
\frac{\phi_{i+1,j,k} - 2\phi_{i,j,k} + \phi_{i-1,j,k}}{\Delta x^2}
+
\frac{\phi_{i,j+1,k} - 2\phi_{i,j,k} + \phi_{i,j-1,k}}{\Delta y^2}
+
\frac{\phi_{i,j,k+1} - 2\phi_{i,j,k} + \phi_{i,j,k-1}}{\Delta z^2}
= 0
$$

Rearranging gives the update used in the iteration:

$$
\phi_{i,j,k}^{\mathrm{new}} =
\frac{
\frac{\phi_W + \phi_E}{\Delta x^2}
+
\frac{\phi_S + \phi_N}{\Delta y^2}
+
\frac{\phi_B + \phi_T}{\Delta z^2}
}{
\frac{2}{\Delta x^2} + \frac{2}{\Delta y^2} + \frac{2}{\Delta z^2}
}
$$

## Iteration strategy

The solver repeatedly updates the field until the maximum absolute change between iterations falls below a tolerance.

That is, it monitors:

$$
\max |\phi^{\mathrm{new}} - \phi|
$$

If this quantity is smaller than the configured tolerance, the iteration stops and the field is treated as converged.

## Output

The example application writes the converged scalar field to a legacy VTK file so it can be visualized in tools such as ParaView.
