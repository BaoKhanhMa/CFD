### Transforming a 3D array into a 1D array

Let's use a standard Cartesian coordinate system (2D first). It would look like this:

        <----------nx---------->
Row 0: [ (0,0) (1,0) (2,0) ... ]
Row 1: [ (0,1) (1,1) (2,1) ... ]
Row 2: [ (0,2) (1,2) (2,2) ... ]

- If we skip one row (increases j), the index increases by nx_
- Therefore, the transformation function: index(i, j) = i + j * nx_

- In 3D, if we increase k, the index increases by nx_ * ny_
- index(i, j, k) = i + j * nx_ k * nx_ * ny_