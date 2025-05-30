# GAULEG.FOR

## Overview

The `GAULEG.FOR` file provides the `GAULEG` subroutine, which calculates the points (abscissas) and weights for Gauss-Legendre quadrature. This numerical integration technique approximates the definite integral of a function by a weighted sum of function values at specific points. The routine calculates these points and weights on a given interval [x1, x2].

This subroutine is functionally almost identical to `GaussLeg` found in `GaussLeg.f90`. It appears to be an older or alternative version, primarily differing in the use of hardcoded Pi, older FORTRAN loop constructs (GOTO), and slight variations in declarations.

## Key Components

*   **SUBROUTINE GAULEG(x1, x2, x, w, n):**
    *   Calculates Gauss-Legendre quadrature abscissas and weights.
    *   **x1 (Input):** REAL*8, the lower bound of the integration interval.
    *   **x2 (Input):** REAL*8, the upper bound of the integration interval.
    *   **x(n) (Output):** REAL*8 array, stores the calculated Gauss-Legendre points (abscissas) scaled to the interval [x1, x2].
    *   **w(n) (Output):** REAL*8 array, stores the calculated Gauss-Legendre weights corresponding to the points in `x(n)`.
    *   **n (Input):** Integer, the number of Gauss-Legendre points and weights to calculate.

## Important Variables/Constants

*   **EPS (Parameter):** DOUBLE PRECISION parameter, set to `3.d-14`. This is the relative precision used as a convergence criterion in Newton's method for finding the roots of the Legendre polynomials.
*   **PI (Hardcoded):** The value `3.141592654d0` is used directly in the formula for the initial guess of the roots.
*   **m:** Integer, `(n+1)/2`. Since the roots of Legendre polynomials are symmetric about 0 on the interval [-1, 1], only half of them need to be explicitly calculated.
*   **xl:** DOUBLE PRECISION, `0.5d0*(x2-x1)`, half the length of the target interval. Used for scaling the roots.
*   **xm:** DOUBLE PRECISION, `0.5d0*(x2+x1)`, the midpoint of the target interval. Used for shifting the roots.
*   **z:** DOUBLE PRECISION, the current approximation of a root of the Legendre polynomial in the interval [-1, 1].
*   **z1:** DOUBLE PRECISION, the previous approximation of `z`, used to check for convergence.
*   **p1, p2, p3:** DOUBLE PRECISION, used in the recurrence relation to evaluate the Legendre polynomial `P_n(z)`. `p1` stores `P_j(z)`, `p2` stores `P_{j-1}(z)`, and `p3` stores `P_{j-2}(z)`.
*   **pp:** DOUBLE PRECISION, the derivative of the Legendre polynomial `P_n(z)`, `P'_n(z)`, used in Newton's method.

## Usage Examples

To integrate a function `F(t)` from `A` to `B` using `N` Gauss-Legendre points:

```fortran
PROGRAM TEST_GAULEG_FOR
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
  CALL GAULEG(A_LIMIT, B_LIMIT, X_GL, W_GL, N_POINTS)

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

END PROGRAM TEST_GAULEG_FOR
```

## Dependencies and Interactions

*   **None:** The `GAULEG` subroutine is self-contained and does not call other external subroutines or functions from this specific project, aside from standard Fortran intrinsics like `ABS`, `COS`.
*   **Numerical Precision:** The use of `REAL*8` / `DOUBLE PRECISION` and a small `EPS` value highlight the need for numerical accuracy.
*   **Symmetry:** The routine exploits the symmetry of Legendre polynomial roots.
*   **Newton's Method:** Uses Newton's method for root-finding.
*   **Comparison with `GaussLeg.f90`:** This routine is very similar to the one in `GaussLeg.f90`. `GaussLeg.f90` uses `ACOS(-1.0D0)` for PI and more modern loop constructs. The choice between them might depend on compiler compatibility or specific project history.
```
