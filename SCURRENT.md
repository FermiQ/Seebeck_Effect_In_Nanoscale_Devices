# SCURRENT.FOR

## Overview

The `SCURRENT.FOR` file provides the `SCURRENT` subroutine, which calculates spin-dependent thermoelectric properties. It computes Seebeck coefficients for spin-up electrons (SEE), spin-down electrons (SEE1), and a combined or spin-difference Seebeck coefficient (SEE2). It also calculates a thermal conductance component (KQ). The method involves integrating spin-up and spin-down transmission functions (interpolated from input data `x1`, `x4up`, `x4dn`) multiplied by Fermi function derivatives and energy terms.

## Key Components

*   **SUBROUTINE SCURRENT(L, TL, TR, DT, x1, x4up, x4dn, SEE, SEE1, SEE2, KQ):**
    *   Calculates spin-dependent Seebeck coefficients and a thermal conductance component.
    *   **L (Input):** Integer, number of data points in the energy grid `x1` and transmission data `x4up`, `x4dn`.
    *   **TL, TR (Input):** Real*8, temperatures of the left and right electron reservoirs (Kelvin).
    *   **DT (Input):** Real*8, (Usage not apparent in this snippet, similar to `CURRENT.FOR`).
    *   **x1(L) (Input):** Real*8 array, energy grid points.
    *   **x4up(L) (Input):** Real*8 array, spin-up transmission probability values.
    *   **x4dn(L) (Input):** Real*8 array, spin-down transmission probability values.
    *   **SEE (Output):** Real*8, Seebeck coefficient for spin-up electrons: `-(CC)/(DD) * 13.60580d0 * 1E6`.
    *   **SEE1 (Output):** Real*8, Seebeck coefficient for spin-down electrons: `-(CC1)/(DD1) * 13.60580D0 * 1E6`.
    *   **SEE2 (Output):** Real*8, Combined/difference Seebeck coefficient: `-(CC-CC1)/(DD+DD1) * 13.60580D0 * 1E6`. (The interpretation might depend on the definition of spin current, whether it's (Up-Down) or (Up+Down) for the denominator).
    *   **KQ (Output):** Real*8, a component of thermal conductance, including contributions from both spins.

## Important Variables/Constants

*   **n (Parameter):** Integer, set to 300. Number of Gauss-Legendre integration points.
*   **COMMON /VL/VL/VR/VR/KT1/KT1/KT2/KT2/VBIAS/VBIAS:** (Same as in `CURRENT.FOR`)
    *   **VL, VR:** Chemical potentials (Rydberg).
    *   **KT1, KT2:** `kB*TL` and `kB*TR` (Rydberg).
*   **y2up(L), y2dn(L):** Real*8 arrays, store second derivatives for spline interpolation of spin-up and spin-down transmissions.
*   **xup, xdn:** Real*8, interpolated transmission values for spin-up and spin-down channels at a given energy `y(j)`.
*   **CC, DD:** Real*8, accumulators for terms related to `K1/T` and `K0` for spin-up electrons, respectively.
*   **CC1, DD1:** Real*8, accumulators for terms related to `K1/T` and `K0` for spin-down electrons, respectively.
*   **K1ULup, K1URup, K2ULup, K2URup:** Real*8, components for spin-up thermal conductance calculation.
*   **K1ULdn, K1URdn, K2ULdn, K2URdn:** Real*8, components for spin-down thermal conductance calculation.
*   **Conversion factors:** `8.6170d-5` (eV/K for Boltzmann const), `13.60580d0` (eV/Rydberg).
*   **Scaling factor:** `151.0d0` used in `K` integral calculations (e.g. `K1ULup`).

## Usage Examples

```fortran
PROGRAM TEST_SCURRENT
  IMPLICIT NONE
  INTEGER, PARAMETER :: L_POINTS = 100
  REAL*8 E_GRID(L_POINTS), T_UP(L_POINTS), T_DOWN(L_POINTS)
  REAL*8 TEMP_L, TEMP_R, DELTA_T
  REAL*8 S_UP, S_DOWN, S_COMB, KQ_OUT
  INTEGER I

  COMMON /VL/VL_C/VR/VR_C/KT1_C/KT1_C/KT2_C/KT2_C/VBIAS_C/VBIAS_C
  REAL*8 VL_C, VR_C, KT1_C, KT2_C, VBIAS_C

! Initialize common block
  VL_C = 0.0D0
  VR_C = 0.0D0
  VBIAS_C = 0.0D0

! Example data
  DO I = 1, L_POINTS
    E_GRID(I) = -0.5D0 + (I-1)*0.01D0
    T_UP(I)   = 0.8D0 * EXP(-(E_GRID(I)-0.05D0)**2 / (2*0.1**2)) ! Spin-up transmission
    T_DOWN(I) = 0.6D0 * EXP(-(E_GRID(I)+0.05D0)**2 / (2*0.1**2)) ! Spin-down transmission
  END DO

  TEMP_L = 300.0D0
  TEMP_R = 301.0D0
  DELTA_T = 0.0D0

  REAL*8 YP1_SPL, YPN_SPL
  YP1_SPL = 1.0D30
  YPN_SPL = 1.0D30

  CALL SCURRENT(L_POINTS, TEMP_L, TEMP_R, DELTA_T, E_GRID, T_UP, T_DOWN, &
                S_UP, S_DOWN, S_COMB, KQ_OUT)

  PRINT *, "Seebeck (Spin Up):", S_UP, "uV/K"
  PRINT *, "Seebeck (Spin Down):", S_DOWN, "uV/K"
  PRINT *, "Seebeck (Combined/Diff):", S_COMB, "uV/K"
  PRINT *, "Thermal Conductance component (KQ):", KQ_OUT

END PROGRAM TEST_SCURRENT
```
*(Note: `FERMIL`, `FERMIR`, `GaussLeg`, `SPLINE`, `SPLINT` functions/subroutines need to be available.)*

## Dependencies and Interactions

*   **GaussLeg (SUBROUTINE):** Called to get Gauss-Legendre integration points and weights. (Note: The first call `CALL GaussLeg(x1(1),x1(L),n,y,w)` integrates over the entire range, while subsequent calls in the loop `CALL GAULEG(x1(i),x1(i+1),y,w,n)` integrate segment-wise. This seems inconsistent. The first call to `GaussLeg` populates `y` and `w` which are then used to calculate `CC, DD, CC1, DD1`. The loop for `K` values re-calls `GAULEG` for each segment.)
*   **SPLINE (SUBROUTINE):** Called to calculate second derivatives for spin-up and spin-down transmission data.
*   **SPLINT (SUBROUTINE):** Called to interpolate transmission values at Gauss-Legendre points.
*   **FERMIL (FUNCTION):** External function, Fermi-Dirac distribution for the left contact.
*   **FERMIR (FUNCTION):** External function, Fermi-Dirac distribution for the right contact.
*   **COMMON /VL/VR/KT1/KT2/VBIAS/:** Shares contact properties.
*   The calculation of `CC, DD, CC1, DD1` uses Gauss-Legendre points/weights derived from the *entire* energy range `x1(1)` to `x1(L)`.
*   The calculation of `K` terms (`K1ULup`, etc.) for `KQ` iterates through segments of the energy range (`x1(i)` to `x1(i+1)`) and calls `GAULEG` for each segment. This is a different integration scheme than for the `CC/DD` terms.
*   The definition of `SEE2` as `-(CC-CC1)/(DD+DD1)` suggests it might be related to a spin polarization of the Seebeck effect or a two-channel model where the conductances (related to DD terms) add.
```
