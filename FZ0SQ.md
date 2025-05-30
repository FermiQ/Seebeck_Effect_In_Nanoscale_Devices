# FZ0SQ.F90

## Overview

The `FZ0SQ.F90` file contains the `FZ0SQ` subroutine. This subroutine calculates a value `F0SQ` (returned as `RESULT`), which appears to be a component of a spectral function or density of states related to surface phonons, specifically for what might be termed the "Zero branch". This routine is called by `SURFPH.F90` as part of its calculation of surface phonon contributions (GG3 term). The calculation involves complex arithmetic and depends on material sound velocities.

## Key Components

*   **SUBROUTINE FZ0SQ(C, RESULT):**
    *   Calculates the `F0SQ` value.
    *   **C (Input):** Double precision. An input parameter, likely a phase velocity or a variable related to it, used in the calculation of `BETA`. In `SURFPH.F90`, `C` takes values from `TC(J)` which are Gaussian quadrature nodes between `VT` and `VL`.
    *   **RESULT (Output):** Double precision. The calculated `F0SQ` value.

## Important Variables/Constants

*   **COMMON/SURF/VL, VT, CR:**
    *   **VL (Input, from common block):** Double precision. Longitudinal sound velocity.
    *   **VT (Input, from common block):** Double precision. Transverse sound velocity.
    *   **CR (Input, from common block):** Double precision. Rayleigh wave velocity.
*   **PI:** Double precision, the mathematical constant Pi, calculated as `ACOS(-1.D0)`.
*   **I:** Complex*16, the imaginary unit `(0.D0, 1.D0)`.
*   **GAMA:** Double precision, `DSQRT(1.D0-(CR/VL)**2)`.
*   **BETA:** Double precision, `DSQRT((C/VT)**2-1.D0)`.
*   **D, E:** Complex*16, intermediate complex variables calculated from `GAMA` and `BETA`. Their expressions are quite involved.
*   **TEMP_P:** Complex*16, `-GAMA*D + I*(1.D0-E)`.
*   **F0SQ:** Double precision, the final calculated value: `DREAL(TEMP_P*DCONJG(TEMP_P))/(2.D0*PI*C*BETA)/C**3`. This represents the squared magnitude of `TEMP_P` scaled by other terms.

## Usage Examples

The `FZ0SQ` subroutine is not typically called directly by the end-user but is a component used by other physics calculation routines. It is called by `SURFPH.F90` within a loop:

```fortran
! Inside SURFPH.F90, during the "ZERO BRANCH" calculation:
! ...
!        DO J=1,N
!                FZ=TC(J)  ! TC(J) are GaussLeg nodes between VT and VL
!                CALL FZ0SQ(FZ,BBB) ! FZ is passed as C to FZ0SQ
!                SP0=SP0+WGTET(J)*BBB
!        END DO
! ...
```
Here, `FZ` (an integration variable related to velocity) is passed as `C` to `FZ0SQ`, and the returned `BBB` (which is `F0SQ`) is used in a weighted sum.

## Dependencies and Interactions

*   **COMMON/SURF/VL, VT, CR:** This common block provides the necessary sound velocities (`VL`, `VT`, `CR`) which must be initialized by a calling routine (e.g., `thermcd.for` which determines `CR`).
*   **SURFPH (SUBROUTINE):** `FZ0SQ` is directly called by `SURFPH.F90` to calculate a component contributing to `GG3`, which is related to the surface phonon spectral function.
*   **Mathematical Functions:** Uses standard Fortran intrinsics like `DSQRT`, `ACOS`, `DREAL`, `DCONJG`.
*   The calculation involves complex number algebra to determine intermediate quantities `D`, `E`, and `TEMP_P`. The final result `F0SQ` is a real number.
```
