# GaussLeg.f90

## Overview

The `GaussLeg.f90` file provides the `GaussLeg` subroutine, which calculates the points (abscissas) and weights for Gauss-Legendre quadrature. This numerical integration technique approximates the definite integral of a function by a weighted sum of function values at specific points. The routine calculates these points and weights on a given interval [x1, x2]. This method is particularly effective for functions that are well-approximated by polynomials. The implementation is based on a routine from "Numerical Recipes in Fortran".

## Key Components

*   **SUBROUTINE GaussLeg(x1, x2, n, x, w):**
    *   Calculates Gauss-Legendre quadrature abscissas and weights.
    *   **x1 (Input):** Real*8, the lower bound of the integration interval.
    *   **x2 (Input):** Real*8, the upper bound of the integration interval.
    *   **n (Input):** Integer, the number of Gauss-Legendre points and weights to calculate.
    *   **x(n) (Output):** Real*8 array, stores the calculated Gauss-Legendre points (abscissas) scaled to the interval [x1, x2].
    *   **w(n) (Output):** Real*8 array, stores the calculated Gauss-Legendre weights corresponding to the points in `x(n)`.

## Important Variables/Constants

*   **EPS (Parameter):** Real*8 parameter, set to `3.d-14`. This is the relative precision used as a convergence criterion in Newton's method for finding the roots of the Legendre polynomials.
*   **pi:** Real*8, the mathematical constant Pi, calculated as `acos(-1.d0)`.
*   **m:** Integer, `(n+1)/2`. Since the roots of Legendre polynomials are symmetric about 0 on the interval [-1, 1], only half of them need to be explicitly calculated.
*   **xl:** Real*8, `0.5*(x2-x1)`, half the length of the target interval. Used for scaling the roots.
*   **xm:** Real*8, `0.5*(x2+x1)`, the midpoint of the target interval. Used for shifting the roots.
*   **z:** Real*8, the current approximation of a root of the Legendre polynomial in the interval [-1, 1].
*   **z1:** Real*8, the previous approximation of `z`, used to check for convergence.
*   **p1, p2, p3:** Real*8, used in the recurrence relation to evaluate the Legendre polynomial `P_n(z)`. `p1` stores `P_j(z)`, `p2` stores `P_{j-1}(z)`, and `p3` stores `P_{j-2}(z)`.
*   **pp:** Real*8, the derivative of the Legendre polynomial `P_n(z)`, `P'_n(z)`, used in Newton's method.

## Usage Examples

To integrate a function `F(t)` from `A` to `B` using `N` Gauss-Legendre points:

```fortran
PROGRAM TEST_GAUSSLEG
  IMPLICIT NONE
  INTEGER, PARAMETER :: N_POINTS = 10
  REAL*8 X_GL(N_POINTS), W_GL(N_POINTS)
  REAL*8 A_LIMIT, B_LIMIT
  REAL*8 INTEGRAL_SUM
  INTEGER I

  EXTERNAL MY_FUNCTION ! The function to be integrated

  A_LIMIT = 0.0D0
  B_LIMIT = 2.0D0

! Get Gauss-Legendre points and weights for the interval [A_LIMIT, B_LIMIT]
  CALL GaussLeg(A_LIMIT, B_LIMIT, N_POINTS, X_GL, W_GL)

  INTEGRAL_SUM = 0.0D0
  DO I = 1, N_POINTS
    INTEGRAL_SUM = INTEGRAL_SUM + W_GL(I) * MY_FUNCTION(X_GL(I))
  END DO

  PRINT *, "Approximate integral of MY_FUNCTION from", A_LIMIT, "to", B_LIMIT
  PRINT *, "using", N_POINTS, "points is:", INTEGRAL_SUM

CONTAINS

  FUNCTION MY_FUNCTION(T_VAL)
    REAL*8 MY_FUNCTION
    REAL*8 T_VAL
    ! Example function: f(t) = t^2
    MY_FUNCTION = T_VAL**2
  END FUNCTION MY_FUNCTION

END PROGRAM TEST_GAUSSLEG
```

## Dependencies and Interactions

*   **None:** The `GaussLeg` subroutine is self-contained and does not call other external subroutines or functions from this specific project, aside from standard Fortran intrinsics like `ABS`, `COS`, `ACOS`.
*   **Numerical Precision:** The comment `High precision is a good idea for this routine` and the use of `Real*8` (double precision) along with a small `EPS` value highlight the need for numerical accuracy in calculating the roots and weights.
*   **Symmetry:** The routine exploits the symmetry of Legendre polynomial roots to reduce computation, calculating roots only for one half of the interval and then mirroring them.
*   **Newton's Method:** The core of the root-finding process for Legendre polynomials uses Newton's method.
*   **Output:** The calculated points `x` and weights `w` are used by other routines that perform numerical integration, such as `SURFPH.F90` in this project.
```
