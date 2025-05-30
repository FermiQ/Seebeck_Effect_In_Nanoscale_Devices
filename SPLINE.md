# SPLINE.FOR

## Overview

The `SPLINE.FOR` file contains the `SPLINE` subroutine, which calculates the second derivatives of a cubic spline interpolating function. Given a set of tabulated data points (x, y), and optionally the first derivatives at the endpoints, this routine computes the second derivatives at each data point. These second derivatives are then used by the `splint` subroutine (from `SPLINT.FOR`) to perform the actual interpolation. The method allows for either specified first derivative boundary conditions or natural spline boundary conditions (zero second derivatives at the ends). This code is adapted from "Numerical Recipes".

## Key Components

*   **SUBROUTINE SPLINE(x, y, n, yp1, ypn, y2):**
    *   Calculates the second derivatives of the interpolating function at the tabulated points.
    *   **x(n) (Input):** A double precision array of x-coordinates of the tabulated data points. Must be in ascending order (`x1 < x2 < ... < xN`).
    *   **y(n) (Input):** A double precision array of y-coordinates of the tabulated data points.
    *   **n (Input):** An integer, the number of data points in `x` and `y`.
    *   **yp1 (Input):** A double precision value for the first derivative at `x(1)`. If `yp1 >= 1.0E30`, a natural spline boundary condition (zero second derivative) is used at `x(1)`.
    *   **ypn (Input):** A double precision value for the first derivative at `x(n)`. If `ypn >= 1.0E30`, a natural spline boundary condition is used at `x(n)`.
    *   **y2(n) (Output):** A double precision array where the calculated second derivatives at each `x(i)` point are stored.

## Important Variables/Constants

*   **NMAX (Parameter):** Integer parameter, set to 2500. It defines the maximum size of the internal temporary array `u`. This implies the subroutine can handle up to 2500 data points.
*   **u(NMAX):** A double precision temporary array used in the algorithm to solve the tridiagonal system for the second derivatives.
*   **sig:** Double precision variable, `(x(i)-x(i-1))/(x(i+1)-x(i-1))`, a ratio of interval lengths.
*   **p:** Double precision variable, `sig*y2(i-1)+2.d0`, an intermediate term in the recurrence relation.
*   **qn, un:** Double precision variables used to handle the boundary condition at `x(n)`.

## Usage Examples

```fortran
PROGRAM TEST_SPLINE
  IMPLICIT NONE
  INTEGER, PARAMETER :: NUM_POINTS = 5
  DOUBLE PRECISION X_VALS(NUM_POINTS), Y_VALS(NUM_POINTS)
  DOUBLE PRECISION Y2_DERIVS(NUM_POINTS)
  DOUBLE PRECISION YP1_BC, YPN_BC
  INTEGER I

! Example data
  DATA X_VALS /1.0D0, 2.0D0, 3.0D0, 4.0D0, 5.0D0/
  DATA Y_VALS /0.5D0, 0.8D0, 0.3D0, 0.6D0, 0.2D0/

! Specify first derivative boundary conditions (e.g., slope of 0 at start, natural at end)
  YP1_BC = 0.0D0
  YPN_BC = 1.0D30 ! Signal for natural spline at the end point

! Calculate second derivatives
  CALL SPLINE(X_VALS, Y_VALS, NUM_POINTS, YP1_BC, YPN_BC, Y2_DERIVS)

  PRINT *, "Calculated second derivatives (y2a):"
  DO I = 1, NUM_POINTS
    PRINT *, "y2(", I, ") = ", Y2_DERIVS(I)
  END DO

! These Y2_DERIVS can now be used with SPLINT for interpolation.
! e.g., CALL splint(X_VALS, Y_VALS, Y2_DERIVS, NUM_POINTS, 2.5D0, Y_INTERP)

END PROGRAM TEST_SPLINE
```

## Dependencies and Interactions

*   **SPLINT (SUBROUTINE):** The output array `y2` from `SPLINE` is essential input for the `splint` subroutine (in `SPLINT.FOR`), which performs the cubic spline interpolation.
*   **Input Data Order:** The input array `x` must be sorted in strictly ascending order.
*   **Boundary Conditions:** The behavior at the endpoints of the spline is controlled by the `yp1` and `ypn` parameters. Large values (>= 0.99d30, effectively 1.0E30) signal the use of natural spline conditions (zero second derivative). Otherwise, the provided values are used as the first derivatives at the boundaries.
*   **NMAX Parameter:** The fixed size `NMAX` for the internal array `u` limits the number of data points that can be processed. If `n` exceeds `NMAX`, an array out-of-bounds error will occur.
```
