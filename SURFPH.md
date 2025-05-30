# SURFPH.F90

## Overview

The `SURFPH.F90` file contains the `SURFPH` subroutine, which calculates components related to surface phonon spectral functions, based on Geller's paper (PRB 64, 155320, 2001). It computes three quantities (GG1, GG2, GG3) representing contributions from different phonon branches (Rayleigh, Plus-Minus, and Zero branches). These calculations are essential for understanding thermal transport at surfaces or in nanostructures where surface effects are significant.

## Key Components

*   **SUBROUTINE SURFPH(GG1, GG2, GG3):**
    *   Calculates three components (GG1, GG2, GG3) of the surface phonon spectral function.
    *   **GG1 (Output):** Double precision. Represents the contribution from the Rayleigh branch, specifically `ALFARSQ`.
    *   **GG2 (Output):** Double precision. Represents the combined contribution from the Plus and Minus branches, scaled by `VL**3*PI`.
    *   **GG3 (Output):** Double precision. Represents the contribution from the Zero branch, scaled by `2.d0*VT**3*PI`.

## Important Variables/Constants

*   **COMMON/SURF/VL, VT, CR:**
    *   **VL:** Double precision input. Longitudinal sound velocity.
    *   **VT:** Double precision input. Transverse sound velocity.
    *   **CR:** Double precision input. Rayleigh wave velocity (presumably calculated elsewhere, e.g., by `thermcd.for` using `ZEROIN`).
*   **N:** Integer, parameter for Gaussian quadrature, set to 2048.
*   **PI:** Double precision, the mathematical constant Pi.
*   **GAMA:** `DSQRT(1.D0-(CR/VL)**2)`, intermediate term for Rayleigh branch calculation.
*   **ETA:** `DSQRT(1.D0-(CR/VT)**2)`, intermediate term for Rayleigh branch calculation.
*   **ALFARSQ:** Calculated value for the Rayleigh branch contribution before scaling `GG1`.
*   **SPRAY:** `ALFARSQ/CR**3`, an intermediate calculation (not directly output).
*   **SC3(2048), WGTES(2048):** Arrays for Gaussian quadrature nodes and weights for Plus-Minus branch calculations.
*   **TC(2048), WGTET(2048):** Arrays for Gaussian quadrature nodes and weights for Zero branch calculations.
*   **SPPLUS, SPMINUS:** Accumulated sums for Plus and Minus branch contributions during numerical integration.
*   **SP0:** Accumulated sum for the Zero branch contribution during numerical integration.

## Usage Examples

The `SURFPH` subroutine is likely called by other routines that require these spectral function components, for example, in the calculation of thermal conductance.

```fortran
PROGRAM TEST_SURFPH
  IMPLICIT DOUBLE PRECISION (A-H,O-Z)
  IMPLICIT INTEGER (I-N)

  COMMON/SURF/VL_COM, VT_COM, CR_COM
  REAL*8 VL_COM, VT_COM, CR_COM
  REAL*8 G1_OUT, G2_OUT, G3_OUT

! Assign realistic values for sound velocities for a material (e.g., Au)
! VL_COM, VT_COM would be material properties.
! CR_COM would typically be solved for using FCR in thermcd.for
  VL_COM = 3.2D5  ! cm/s
  VT_COM = 1.2D5  ! cm/s
! Assume CR has been found, e.g. CR = 0.9 * VT_COM
  CR_COM = 0.9D0 * VT_COM

  CALL SURFPH(G1_OUT, G2_OUT, G3_OUT)

  PRINT *, "GG1 (Rayleigh Contribution Factor):", G1_OUT
  PRINT *, "GG2 (Plus-Minus Branch Factor):", G2_OUT
  PRINT *, "GG3 (Zero Branch Factor):", G3_OUT

END PROGRAM TEST_SURFPH
```

## Dependencies and Interactions

*   **COMMON/SURF/VL, VT, CR:** This common block is used to receive the longitudinal sound velocity (VL), transverse sound velocity (VT), and Rayleigh wave velocity (CR) as inputs. These values are typically defined and potentially calculated (for CR) in a calling routine like `thermcd`.
*   **GaussLeg (SUBROUTINE):** This external subroutine is called to generate Gaussian quadrature nodes (stored in `SC3` and `TC`) and weights (stored in `WGTES` and `WGTET`). This is used for numerical integration of terms related to the Plus-Minus and Zero branches.
*   **FZPNSQ (SUBROUTINE):** This external subroutine is called within the loop for the Plus-Minus branch calculation. It takes `F1Z` (which is `1.0/SC3(J)`) as input and returns `FPSQ` and `FNSQ`.
*   **FZ0SQ (SUBROUTINE):** This external subroutine is called within the loop for the Zero branch calculation. It takes `FZ` (which is `TC(J)`) as input and returns `BBB`.
*   The results `GG1`, `GG2`, and `GG3` are passed back to the calling routine (e.g., `thermcd`) and are used in subsequent calculations, such as determining thermal capacitance (`Cth` in `thermcd.for`).
```
