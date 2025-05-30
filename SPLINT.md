# SPLINT.FOR

## Overview

The `SPLINT.FOR` file provides a subroutine `splint` that performs cubic spline interpolation. Given a set of tabulated data points (x, y) and their second derivatives (calculated by a routine like `spline`), this subroutine computes the interpolated y-value for a given x-value. This is a standard numerical technique used for interpolating data when a smooth curve representation is needed. The code is adapted from "Numerical Recipes in FORTRAN 77".

## Key Components

*   **SUBROUTINE splint(xa, ya, y2a, n, x, y):**
    *   Performs cubic spline interpolation.
    *   **xa(n) (Input):** A real*8 array of x-coordinates of the tabulated data points. Must be in ascending order.
    *   **ya(n) (Input):** A real*8 array of y-coordinates of the tabulated data points.
    *   **y2a(n) (Input):** A real*8 array containing the second derivatives of the interpolating function at the `xa` points. This is typically pre-calculated by a `spline` subroutine.
    *   **n (Input):** An integer, the number of data points in `xa`, `ya`, and `y2a`.
    *   **x (Input):** A real*8 value, the x-coordinate at which to interpolate.
    *   **y (Output):** A real*8 value, the interpolated y-coordinate corresponding to `x`.

## Important Variables/Constants

*   **klo, khi:** Integers used as lower and upper indices for binary searching the interval in `xa` that brackets the input `x`.
*   **k:** Integer, the midpoint index during the binary search.
*   **h:** Real*8, the width of the bracketing interval (`xa(khi) - xa(klo)`).
*   **a, b:** Real*8, weighting factors used in the cubic spline interpolation formula, calculated as `(xa(khi)-x)/h` and `(x-xa(klo))/h` respectively.

## Usage Examples

To use `splint`, you first need to calculate the second derivatives `y2a` using a `spline` routine (which is defined in `SPLINE.FOR` in this project).

```fortran
PROGRAM TEST_SPLINT
  IMPLICIT NONE
  INTEGER, PARAMETER :: NUM_POINTS = 5
  REAL*8 X_VALS(NUM_POINTS), Y_VALS(NUM_POINTS), Y2_DERIVS(NUM_POINTS)
  REAL*8 X_INTERP, Y_INTERP
  INTEGER I

! Example data
  DATA X_VALS /1.0, 2.0, 3.0, 4.0, 5.0/
  DATA Y_VALS /0.5, 0.8, 0.3, 0.6, 0.2/

! Placeholder for y2_derivs - these would be computed by SPLINE subroutine
! For example: CALL spline(X_VALS, Y_VALS, NUM_POINTS, 1.0E30, 1.0E30, Y2_DERIVS)
! Using arbitrary values for demonstration here if SPLINE is not directly callable
  DATA Y2_DERIVS /0.0, -0.5, 0.4, -0.3, 0.0/ ! Replace with actual spline output

! Point at which to interpolate
  X_INTERP = 2.5

! Perform spline interpolation
  CALL splint(X_VALS, Y_VALS, Y2_DERIVS, NUM_POINTS, X_INTERP, Y_INTERP)

  PRINT *, "Interpolated value at x =", X_INTERP, " is y =", Y_INTERP

END PROGRAM TEST_SPLINT
```

## Dependencies and Interactions

*   **SPLINE (SUBROUTINE):** The `splint` subroutine relies on pre-calculated second derivative values (`y2a`). These values are typically computed by a companion subroutine named `spline` (found in `SPLINE.FOR` in this repository), which takes the `xa` and `ya` arrays as input.
*   **Input Data Order:** The input array `xa` must be sorted in ascending order for the binary search algorithm to function correctly.
*   **Error Handling:** The subroutine includes a basic error check: `if (h.eq.0.) pause 'bad xa input in splint'`. This occurs if two consecutive `xa` values are identical, which would lead to division by zero.
```
