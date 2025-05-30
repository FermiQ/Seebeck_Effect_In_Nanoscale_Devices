# ZEROIN.F90

## Overview

The `ZEROIN.F90` file contains the `ZEROIN` subroutine, which is designed to find a zero of a continuous function within a given interval [A, B]. It requires that the function values at the endpoints, FCT(A) and FCT(B), have opposite signs. The method is a combination of bisection and secant methods, incorporating inverse quadratic interpolation to efficiently locate the zero. This subroutine is a numerical method crucial for root finding in various computational contexts within the project.

## Key Components

*   **SUBROUTINE ZEROIN(FCT, ABSERR, RELERR, MAXFCT, ITAPE, A, B, FB, NUMFCT, IERR):**
    *   **FCT:** The input function (must be `EXTERNAL` or `INTRINSIC`) for which the zero is to be found.
    *   **ABSERR, RELERR:** Absolute and relative error tolerances for the zero.
    *   **MAXFCT:** Maximum number of function evaluations allowed.
    *   **ITAPE:** Unit number for an optional output data set for intermediate results (if > 0).
    *   **A, B:** Input interval endpoints. On output, B is the approximate zero.
    *   **FB:** Output function value at the approximate zero B.
    *   **NUMFCT:** Output number of function evaluations performed.
    *   **IERR:** Output error parameter indicating the status of the zero-finding process.

## Important Variables/Constants

*   **EPS:** Machine epsilon (calculated as four times the machine constant `FMACHP`). Used to adjust input error bounds `ABSERR` and `RELERR` if they are too small.
*   **FMACHP:** Machine precision constant, dynamically calculated.
*   **FA, FB, FC:** Function values at points A, B, and C respectively. C is an auxiliary point used in the algorithm.
*   **XM:** Half the interval length (C - B)/2, used in convergence checks.
*   **TOL1:** Tolerance value calculated from `ABSERR` and `RELERR` for the convergence criterion: `0.5D0 * (ABSERR + RELERR * DABS(B))`.
*   **IERR values:**
    *   -2: Invalid input parameters (e.g., `ABSERR` or `RELERR` negative, `MAXFCT < 1`).
    *   -1: `FCT(A)*FCT(B) >= 0`, meaning no sign change detected, so a zero might not be bracketed.
    *   0: `A` or `B` is a numerical zero of `FCT`.
    *   1: `B` is an exact zero with `FCT(B)=0.0`.
    *   2: Desired accuracy has been achieved (`DABS(XM) <= TOL1`).
    *   3: Maximum number of function evaluations (`MAXFCT`) reached without achieving desired accuracy.

## Usage Examples

```fortran
! Assuming FCT is a defined DOUBLE PRECISION FUNCTION FCT(X)
EXTERNAL FCT
DOUBLE PRECISION A_VAL, B_VAL, ABS_ERR, REL_ERR, FUNC_B, RESULT_B
INTEGER MAX_EVAL, TAPE_NO, NUM_EVAL, ERR_CODE

A_VAL = 1.0D0
B_VAL = 2.0D0
ABS_ERR = 1.0D-8
REL_ERR = 1.0D-8
MAX_EVAL = 100
TAPE_NO = 0 ! No intermediate output

CALL ZEROIN(FCT, ABS_ERR, REL_ERR, MAX_EVAL, TAPE_NO, A_VAL, RESULT_B, FUNC_B, NUM_EVAL, ERR_CODE)

IF (ERR_CODE .GT. 0) THEN
  PRINT *, "Approximate zero found at:", RESULT_B
  PRINT *, "Function value at zero:", FUNC_B
  PRINT *, "Number of evaluations:", NUM_EVAL
ELSE
  PRINT *, "Error in ZEROIN, code:", ERR_CODE
ENDIF
```

## Dependencies and Interactions

*   **MACHPD (INTEGER FUNCTION):** This function is called to determine machine precision. It's likely an external or system-provided function. The comment `Required subroutines: MACHPD, BI` indicates this dependency.
*   **BI (SUBROUTINE):** This subroutine is called under certain conditions (e.g., when interpolation steps are not effective or when the interval needs to be bisected). The comment `Required subroutines: MACHPD, BI` indicates this dependency.
*   **FCT (EXTERNAL DOUBLE PRECISION FUNCTION):** The user must provide this function, representing the equation for which the root is sought.
*   The subroutine may write intermediate results to a data set specified by `ITAPE` if `ITAPE > 0`.
```
