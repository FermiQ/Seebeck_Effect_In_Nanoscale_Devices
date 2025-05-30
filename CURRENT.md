# CURRENT.FOR

## Overview

The `CURRENT.FOR` file defines the `CURRENT` subroutine, which calculates various thermoelectric properties for a system, likely a molecular junction. These properties include electrical current (EG), heat current (HG), Seebeck coefficient (SEE), and thermal conductance (KQ) components. The calculation is performed by integrating the transmission function (interpolated from input data `x1`, `x4`) multiplied by Fermi function differences and energy-dependent terms over an energy range. This version seems to be for a spinless case.

## Key Components

*   **SUBROUTINE CURRENT(L, TL, TR, DT, x1, x4, CC, DD, EG, HG, HH, SEE, KQ):**
    *   Calculates electrical and heat currents, Seebeck coefficient, and thermal conductance terms.
    *   **L (Input):** Integer, the number of data points in the input energy grid `x1` and transmission data `x4`.
    *   **TL, TR (Input):** Real*8, temperatures of the left and right electron reservoirs (in Kelvin).
    *   **DT (Input):** Real*8, likely a time step or temperature difference, though its direct usage in the main calculations isn't apparent in this snippet (it's checked in a `do while (DT.LT.10.0d0)` loop that is commented out).
    *   **x1(L) (Input):** Real*8 array, energy grid points.
    *   **x4(L) (Input):** Real*8 array, transmission probability values corresponding to the energy grid `x1`.
    *   **CC (Output):** Real*8, accumulates terms related to `K1/T` for Seebeck calculation. `K1 = - integral ( (E-Ef)*f*(1-f)/(kB*T) * Transmission(E) dE )`.
    *   **DD (Output):** Real*8, accumulates terms related to `K0` for Seebeck calculation. `K0 = - integral ( f*(1-f)/(kB*T) * Transmission(E) dE )`.
    *   **EG (Output):** Real*8, calculated electrical current, scaled by `1/151.0d0`. `I = integral ( Transmission(E) * (f_L(E) - f_R(E)) dE )`.
    *   **HG (Output):** Real*8, calculated heat current. `J_Q = integral ( Transmission(E) * (E - mu_L) * (f_L(E) - f_R(E)) dE )`. (Note: `mu_L` is `VL` here).
    *   **HH (Output):** Real*8, accumulates terms similar to `DD` but scaled by `1/151.0d0`. Related to thermal averaging.
    *   **SEE (Output):** Real*8, calculated Seebeck coefficient `S = -(CC)/(DD) * 13.60580d0 * 1E6` (in microVolts/K).
    *   **KQ (Output):** Real*8, a component of thermal conductance.

## Important Variables/Constants

*   **n (Parameter):** Integer, set to 200. Number of Gauss-Legendre integration points used between each pair of input energy grid points.
*   **COMMON /VL/VL/VR/VR/KT1/KT1/KT2/KT2/VBIAS/VBIAS:**
    *   **VL, VR:** Real*8, chemical potentials (or effective Fermi levels) of the left and right contacts.
    *   **KT1, KT2:** Real*8, `kB*TL` and `kB*TR` respectively, in Rydberg units (where kB is Boltzmann constant).
    *   **VBIAS:** Real*8, bias voltage (not explicitly used in the calculations shown in this snippet but part of the common block).
*   **y2(L):** Real*8 array, stores second derivatives for spline interpolation, calculated from `x1` and `x4`.
*   **y(n), w(n):** Real*8 arrays, store abscissas and weights for Gauss-Legendre integration.
*   **AA:** Real*8, `FERMIR(y(j)) - FERMIL(y(j))`, difference of Fermi functions at energy `y(j)`.
*   **K1UL, K1UR, K2UL, K2UR:** Real*8, components used in calculating `KQ`, related to integrals involving powers of `(E - mu)` and Fermi function derivatives.
    *   `K1UL/R ~ integral (E-mu_L/R) * Trans(E) * df/dE_L/R dE`
    *   `K2UL/R ~ integral (E-mu_L/R)^2 * Trans(E) * df/dE_L/R dE`
*   **Conversion factors:** `8.6170d-5` (eV/K for Boltzmann const), `13.60580d0` (eV/Rydberg).
*   **Scaling factor:** `151.0d0` appears in scaling `EG` and terms for `K1UL, K1UR, K2UL, K2UR`. Its physical meaning is unclear from the code alone.

## Usage Examples

```fortran
PROGRAM SIMULATE_CURRENT
  IMPLICIT NONE
  INTEGER, PARAMETER :: L_POINTS = 100
  REAL*8 ENERGY_GRID(L_POINTS), TRANSMISSION(L_POINTS)
  REAL*8 TEMP_L, TEMP_R, DELTA_T
  REAL*8 CC_OUT, DD_OUT, EG_OUT, HG_OUT, HH_OUT, SEE_OUT, KQ_OUT
  INTEGER I

  COMMON /VL/VL_COM/VR/VR_COM/KT1_COM/KT1_COM/KT2_COM/KT2_COM/VBIAS_COM/VBIAS_COM
  REAL*8 VL_COM, VR_COM, KT1_COM, KT2_COM, VBIAS_COM

! Initialize common block variables
  VL_COM = 0.0D0  ! Fermi level of left contact (Rydberg)
  VR_COM = 0.0D0  ! Fermi level of right contact (Rydberg) - for zero bias
  VBIAS_COM = 0.0D0

! Initialize energy grid and transmission data (example)
  DO I = 1, L_POINTS
    ENERGY_GRID(I) = -0.5D0 + (I-1)*0.01D0 ! Example energy range
    TRANSMISSION(I) = EXP(-(ENERGY_GRID(I))**2 / (2*0.1**2)) ! Example Gaussian transmission
  END DO

  TEMP_L = 300.0D0 ! Kelvin
  TEMP_R = 300.0D0 ! Kelvin
  DELTA_T = 0.0D0  ! Not actively used here

! yp1, ypn for spline (typically large for natural spline, or specific derivatives)
  REAL*8 YP1_SPLINE, YPN_SPLINE
  YP1_SPLINE = 1.0D30
  YPN_SPLINE = 1.0D30


  CALL CURRENT(L_POINTS, TEMP_L, TEMP_R, DELTA_T, ENERGY_GRID, TRANSMISSION, &
               CC_OUT, DD_OUT, EG_OUT, HG_OUT, HH_OUT, SEE_OUT, KQ_OUT)

  PRINT *, "Calculated Current (EG):", EG_OUT
  PRINT *, "Calculated Heat Current (HG):", HG_OUT
  PRINT *, "Seebeck Coefficient (SEE):", SEE_OUT, "uV/K"
  PRINT *, "Thermal Conductance component (KQ):", KQ_OUT
  PRINT *, "CC, DD, HH:", CC_OUT, DD_OUT, HH_OUT

END PROGRAM SIMULATE_CURRENT
```
*(Note: `FERMIL` and `FERMIR` functions need to be available for this example)*

## Dependencies and Interactions

*   **GAULEG (SUBROUTINE):** Called to get Gauss-Legendre integration points and weights.
*   **SPLINE (SUBROUTINE):** Called to calculate second derivatives of the transmission data for interpolation.
*   **SPLINT (SUBROUTINE):** Called to interpolate transmission values at Gauss-Legendre points.
*   **FERMIL (FUNCTION):** External function, calculates the Fermi-Dirac distribution `f(E, mu_L, T_L)`.
*   **FERMIR (FUNCTION):** External function, calculates the Fermi-Dirac distribution `f(E, mu_R, T_R)`.
*   **COMMON /VL/VR/KT1/KT2/VBIAS/:** Shares contact chemical potentials, `kB*T` values, and bias voltage with other parts of the program (e.g., `FERMIL`, `FERMIR`).
*   The input arrays `x1` (energy) and `x4` (transmission) define the system's basic transport characteristics.
*   The subroutine performs numerical integration over energy, breaking the total energy range (defined by `x1(1)` to `x1(L)`) into `L-1` segments, and applying `n`-point Gauss-Legendre quadrature to each segment.
```
