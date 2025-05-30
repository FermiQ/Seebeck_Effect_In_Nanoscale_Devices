# CURRENT_EP.FOR

## Overview

The `CURRENT_EP.FOR` file defines the `CURRENT_EP` subroutine, which calculates thermoelectric properties for a system, likely a molecular junction, while explicitly considering electron-phonon (EP) interactions. Similar to the `CURRENT` subroutine, it computes electrical current (EG), heat current (HG), and related Seebeck coefficients (SEE, SEE1). However, it includes terms (`CC1F`, `DD1F`) that account for inelastic scattering processes involving phonon absorption and emission, summed over different phonon modes. The calculation integrates the transmission function (interpolated from `x1`, `x4`) multiplied by Fermi functions and energy-dependent terms, now also including phonon occupation numbers (`AVG(k)`) and coupling strengths (`CLR`, `CRR`).

## Key Components

*   **SUBROUTINE CURRENT_EP(L, TL, TR, DT, fr, x1, x4, CC, DD, EG, HG, HH, SEE, KQ, CLR, CRR, Tw, SEE1):**
    *   Calculates electrical and heat currents, Seebeck coefficients, and thermal conductance terms, including electron-phonon interaction effects.
    *   **L (Input):** Integer, number of data points in the energy grid `x1` and transmission data `x4`.
    *   **TL, TR (Input):** Real*8, temperatures of the left and right electron reservoirs (Kelvin).
    *   **DT (Input):** Real*8, (Usage not apparent in this snippet, similar to `CURRENT.FOR`).
    *   **fr(NMODE) (Input):** Real*8 array, frequencies of the phonon modes (in Rydberg units, as `AKB` is later used with it).
    *   **x1(L) (Input):** Real*8 array, energy grid points.
    *   **x4(L) (Input):** Real*8 array, transmission probability values (elastic part) corresponding to `x1`.
    *   **CC (Output):** Real*8, elastic contribution to `K1/T` terms for Seebeck.
    *   **DD (Output):** Real*8, elastic contribution to `K0` terms for Seebeck.
    *   **EG (Output):** Real*8, electrical current (elastic part, scaled by `1/151.0d0`).
    *   **HG (Output):** Real*8, heat current (elastic part).
    *   **HH (Output):** Real*8, elastic contribution to thermal averaging term, scaled by `1/151.0d0`.
    *   **SEE (Output):** Real*8, Seebeck coefficient based on elastic contributions: `-(CC)/(DD) * 13.60580d0 * 1E6`.
    *   **KQ (Output):** Real*8, a component of thermal conductance, calculated from `K0LL, K0LR, K1UL, K1UR, K2UL, K2UR` which seem to be elastic contributions.
    *   **CLR(NMODE), CRR(NMODE) (Input):** Real*8 arrays, likely representing coupling strengths of electrons to phonon mode `k` for left and right contacts/processes.
    *   **Tw (Input):** Real*8, temperature of the phonon bath (local temperature of the junction/molecule, in Kelvin).
    *   **SEE1 (Output):** Real*8, Seebeck coefficient including electron-phonon effects: `-(CC-CC1F)/(DD-DD1F)*13.6058d0*1E6`.
    *   **CC1F (Implicit Output via accumulation):** Real*8, sum over modes of `CC1`, the inelastic (phonon-mediated) contribution to `K1/T`.
    *   **DD1F (Implicit Output via accumulation):** Real*8, sum over modes of `DD1`, the inelastic (phonon-mediated) contribution to `K0`.

## Important Variables/Constants

*   **n (Parameter):** Integer, set to 20. Number of Gauss-Legendre integration points per segment.
*   **NMODE (Parameter):** Integer, set to 3. Number of phonon modes considered.
*   **AKB (Parameter):** Real*8, Boltzmann constant `8.617d-5` (eV/K). Note: `fr(k)/(AKB*Tw)` suggests `fr(k)` should be in eV for this expression, or `AKB` should be in Rydberg/K if `fr(k)` is in Rydberg. Given other conversions, `fr(k)` is likely in Rydberg.
*   **COMMON /VL/VR/KT1/KT1/KT2/KT2/VBIAS/VBIAS:** (Same as in `CURRENT.FOR`)
    *   **VL, VR:** Chemical potentials (Rydberg).
    *   **KT1, KT2:** `kB*TL` and `kB*TR` (Rydberg).
*   **AVG(NMODE):** Real*8 array, Bose-Einstein distribution `1/(exp(hbar*omega_k / kB*Tw) - 1)` for each phonon mode `k` at temperature `Tw`.
*   **CC1, DD1:** Real*8, temporary accumulators for inelastic contributions for a single mode before being added to `CC1F`, `DD1F`. These expressions are complex, involving products of Fermi functions at energies `y(j)`, `y(j)+fr(k)`, `y(j)-fr(k)`, and coupling terms `CLRc`, `CRRc`.
*   **K0LL, K0LR, K1UL, K1UR, K2UL, K2UR:** (Similar to `CURRENT.FOR`) Components for calculating `KQ`, seemingly based on elastic scattering.
*   **TEMPL, TEMPR:** `kB*TL` and `kB*TR` in Rydberg, used within the integration loops.
*   **Conversion factors:** `13.6058` (eV/Rydberg).
*   **Scaling factor:** `151.0d0` (Same as in `CURRENT.FOR`).

## Usage Examples

```fortran
PROGRAM SIMULATE_CURRENT_EP
  IMPLICIT NONE
  INTEGER, PARAMETER :: L_PTS = 100, NUM_MODES = 3
  REAL*8 E_GRID(L_PTS), T_ELASTIC(L_PTS)
  REAL*8 T_LEFT, T_RIGHT, T_PHONON, DELTA_TEMP
  REAL*8 PHONON_FREQ(NUM_MODES), COUPL_L(NUM_MODES), COUPL_R(NUM_MODES)
  REAL*8 CC_OUT, DD_OUT, EG_OUT, HG_OUT, HH_OUT, SEE_ELAST_OUT, KQ_OUT, SEE_INELAST_OUT

  COMMON /VL/VL_C/VR/VR_C/KT1_C/KT1_C/KT2_C/KT2_C/VBIAS_C/VBIAS_C
  REAL*8 VL_C, VR_C, KT1_C, KT2_C, VBIAS_C

! Initialize common block
  VL_C = 0.0D0
  VR_C = 0.0D0
  VBIAS_C = 0.0D0

! Example data
  INTEGER I
  DO I = 1, L_PTS
    E_GRID(I) = -0.5D0 + (I-1)*0.01D0
    T_ELASTIC(I) = EXP(-(E_GRID(I))**2 / (2*0.1**2))
  END DO

  T_LEFT = 300.0D0
  T_RIGHT = 301.0D0
  T_PHONON = 300.5D0
  DELTA_TEMP = 0.0D0 ! Not actively used

! Phonon properties (example - ensure units are consistent, likely Rydberg for frequencies)
  PHONON_FREQ(1) = 0.01D0 / 13.6058 ! approx 0.01 eV
  PHONON_FREQ(2) = 0.02D0 / 13.6058 ! approx 0.02 eV
  PHONON_FREQ(3) = 0.03D0 / 13.6058 ! approx 0.03 eV
  DO I = 1, NUM_MODES
    COUPL_L(I) = 0.01D0 ! Example coupling strength
    COUPL_R(I) = 0.01D0 ! Example coupling strength
  END DO

  REAL*8 YP1_SPL, YPN_SPL
  YP1_SPL = 1.0D30
  YPN_SPL = 1.0D30

  CALL CURRENT_EP(L_PTS, T_LEFT, T_RIGHT, DELTA_TEMP, PHONON_FREQ, E_GRID, T_ELASTIC, &
                  CC_OUT, DD_OUT, EG_OUT, HG_OUT, HH_OUT, SEE_ELAST_OUT, KQ_OUT, &
                  COUPL_L, COUPL_R, T_PHONON, SEE_INELAST_OUT)

  PRINT *, "Elastic Seebeck (SEE):", SEE_ELAST_OUT, "uV/K"
  PRINT *, "Inelastic Seebeck (SEE1):", SEE_INELAST_OUT, "uV/K"
  PRINT *, "Electrical Current (EG):", EG_OUT
  PRINT *, "Heat Current (HG):", HG_OUT
  PRINT *, "KQ:", KQ_OUT

END PROGRAM SIMULATE_CURRENT_EP
```
*(Note: `FERMIL` and `FERMIR` functions need to be available.)*

## Dependencies and Interactions

*   **GAULEG (SUBROUTINE):** Called for Gauss-Legendre integration points/weights.
*   **SPLINE (SUBROUTINE):** Called to calculate second derivatives for interpolating elastic transmission.
*   **SPLINT (SUBROUTINE):** Called to interpolate elastic transmission values.
*   **FERMIL (FUNCTION):** External function, Fermi-Dirac distribution for the left contact.
*   **FERMIR (FUNCTION):** External function, Fermi-Dirac distribution for the right contact.
*   **COMMON /VL/VR/KT1/KT2/VBIAS/:** Shares contact properties.
*   The subroutine iterates `NMODE` times, once for each phonon mode, accumulating the inelastic contributions (`CC1F`, `DD1F`) to the Seebeck coefficient.
*   The elastic part of the calculation (`EG`, `HG`, `CC`, `DD`, `HH`, `KQ`) is recalculated inside the loop over phonon modes. This seems redundant for `EG`, `HG`, `HH`, `CC`, `DD` if `x1,x4,TL,TR` don't change with mode `k`. However, `CC` and `DD` are re-initialized for each mode and then used for the elastic `SEE` calculation after the loop, which means only the last mode's elastic `CC,DD` are used for `SEE` if not careful. This might be a bug or a misunderstanding of the intended logic for `SEE`. The final `SEE` is calculated using the `CC` and `DD` from the *last* phonon mode iteration's elastic calculation. `SEE1` uses the accumulated `CC1F` and `DD1F` and the final `CC` and `DD`.
```
