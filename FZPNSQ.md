# FZPNSQ.F90

## Overview

The `FZPNSQ.F90` file contains the `FZPNSQ` subroutine. This subroutine calculates two values, `FPSQ` and `FNSQ`. These quantities are likely components of a surface phonon spectral function, corresponding to what might be termed "Plus" and "Minus" branches (or P and N waves). This routine is called by `SURFPH.F90` as part of its calculation of surface phonon contributions (specifically, the GG2 term). The calculations involve complex arithmetic and depend on material sound velocities.

## Key Components

*   **SUBROUTINE FZPNSQ(C, FPSQ, FNSQ):**
    *   Calculates the `FPSQ` and `FNSQ` values.
    *   **C (Input):** Double precision. An input parameter, likely related to a wave vector or inverse phase velocity. In `SURFPH.F90`, `C` is `1.0/SC3(J)`, where `SC3(J)` are Gaussian quadrature nodes for integrating over wavevector-like variables.
    *   **FPSQ (Output):** Double precision. The calculated value for the "Plus" branch component.
    *   **FNSQ (Output):** Double precision. The calculated value for the "Minus" branch component.

## Important Variables/Constants

*   **COMMON/SURF/VL, VT, CR:**
    *   **VL (Input, from common block):** Double precision. Longitudinal sound velocity.
    *   **VT (Input, from common block):** Double precision. Transverse sound velocity.
    *   **CR (Input, from common block):** Double precision. Rayleigh wave velocity. (Although `CR` is in the common block, it is not explicitly used in the calculations within `FZPNSQ`).
*   **PI:** Double precision, the mathematical constant Pi, calculated as `ACOS(-1.D0)`.
*   **I:** Complex*16, the imaginary unit `(0.D0, 1.D0)`.
*   **ALFA:** Double precision, `DSQRT((C/VL)**2-1.D0)`. Note: If `C` is an inverse velocity `1/v_phase`, then `C/VL` would be `v_sound_L / v_phase`. This term requires `C > VL` for `ALFA` to be real, or `(C/VL)**2 > 1`. Given the context in `SURFPH.F90` where `C` is `1/SC3(J)` and `SC3(J)` are integration points from `0` to `1/VL`, `C` would range from `VL` to infinity.
*   **BETA:** Double precision, `DSQRT((C/VT)**2-1.D0)`. Similar to `ALFA`, this requires `C > VT`.
*   **A, B:** Double precision, intermediate variables calculated using `ALFA` and `BETA`.
    *   `A=((BETA**2-1.D0)**2-4.D0*ALFA*BETA) / ((BETA**2-1.D0)**2+4.D0*ALFA*BETA)`
    *   `B=(4.D0*DSQRT(ALFA*BETA)*(BETA**2-1.D0)) / ((BETA**2-1.D0)**2+4.D0*ALFA*BETA)`
*   **TEMP_P:** Complex*16, `DSQRT(ALFA)*(1.D0+A+I*B)+I*(1.D0-A-I*B)/DSQRT(BETA)`.
*   **TEMP_N:** Complex*16, `-DSQRT(ALFA)*(1.D0+A-I*B)+I*(1.D0-A+I*B)/DSQRT(BETA)`.
*   **FPSQ (Final Value):** `DREAL(TEMP_P*DCONJG(TEMP_P))/(4.D0*PI*C)/C**3`. This represents the squared magnitude of `TEMP_P` scaled by other terms.
*   **FNSQ (Final Value):** `DREAL(TEMP_N*DCONJG(TEMP_N))/(4.D0*PI*C)/C**3`. This represents the squared magnitude of `TEMP_N` scaled by other terms.

## Usage Examples

The `FZPNSQ` subroutine is an internal component used by other physics calculation routines, specifically `SURFPH.F90`.

```fortran
! Inside SURFPH.F90, during the "PLUS-MINUS BRANCH" calculation:
! ...
!        DO J=1,N
!                F1Z=1.D0/SC3(J) ! SC3(J) are GaussLeg nodes from X1=0 to X2=1/VL
!                CALL FZPNSQ(F1Z,FPSQ_VAL,FNSQ_VAL) ! F1Z is passed as C
!                SPPLUS=SPPLUS+(FPSQ_VAL/(SC3(J)**2))*WGTES(J)
!                SPMINUS=SPMINUS+(FNSQ_VAL/(SC3(J)**2))*WGTES(J)
!        END DO
! ...
```
Here, `F1Z` (which has units of velocity if `SC3(J)` has units of inverse velocity) is passed as `C` to `FZPNSQ`. The returned `FPSQ_VAL` and `FNSQ_VAL` are then used in weighted sums.

## Dependencies and Interactions

*   **COMMON/SURF/VL, VT, CR:** This common block provides `VL` and `VT`. These sound velocities must be initialized by a calling routine (e.g., `thermcd.for`). `CR` is present but not used.
*   **SURFPH (SUBROUTINE):** `FZPNSQ` is directly called by `SURFPH.F90` to calculate components contributing to `GG2`, related to the surface phonon spectral function.
*   **Mathematical Functions:** Uses standard Fortran intrinsics like `DSQRT`, `ACOS`, `DREAL`, `DCONJG`.
*   The calculations for `ALFA` and `BETA` involve `(C/VL)**2-1.D0` and `(C/VT)**2-1.D0`. If `C` can be less than `VL` or `VT`, these terms can become negative, leading to errors with `DSQRT` if not handled (e.g. by ensuring `C` is always greater or by using complex square roots if the physics implies evanescent modes). Given the integration limits in `SURFPH` for `SC3(J)` (from 0 to `1/VL`), `C = 1/SC3(J)` will range from `VL` to infinity, so `(C/VL)**2 >= 1` and `(C/VT)**2 >= (VL/VT)**2`. If `VL > VT`, then `(C/VT)**2` will also be `>=1`.
*   The routine involves complex number algebra for `TEMP_P` and `TEMP_N`. The final results `FPSQ` and `FNSQ` are real numbers.
```
