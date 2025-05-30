# FERMIR.FOR

## Overview

The `FERMIR.FOR` file defines the `FERMIR` function, which calculates the Fermi-Dirac distribution function for the right electron reservoir/contact. The Fermi-Dirac distribution gives the probability that an electron state at a given energy `E` is occupied, based on the chemical potential of the right contact (`VR`) and its thermal energy (`KT2 = kB*TR`).

## Key Components

*   **REAL*8 FUNCTION FERMIR(E):**
    *   Calculates the Fermi-Dirac distribution `f_R(E) = 1 / (1 + exp((E - VR) / KT2))`.
    *   **E (Input):** Real*8, the energy at which to evaluate the distribution.
    *   **Returns (FERMIR):** Real*8, the value of the Fermi-Dirac distribution function.

## Important Variables/Constants

*   **COMMON /VL/VL/VR/VR/KT1/KT1/KT2/KT2:**
    *   **VR (Input, from common block):** Real*8, the chemical potential (Fermi level) of the right contact.
    *   **KT2 (Input, from common block):** Real*8, the thermal energy of the right contact (`kB * TR`, where `kB` is the Boltzmann constant and `TR` is the temperature of the right contact). This value is expected to be in the same energy units as `E` and `VR`.

## Usage Examples

The `FERMIR` function is typically called within integration loops in routines like `CURRENT.FOR`, `CURRENT_EP.FOR`, and `SCURRENT.FOR` to determine the occupation probability of states at different energies in the right contact.

```fortran
PROGRAM TEST_FERMIR
  IMPLICIT NONE
  REAL*8 ENERGY_VAL, FERMI_DIST_R

  COMMON /VL/VL_C/VR/VR_C/KT1_C/KT1_C/KT2_C/KT2_C
  REAL*8 VL_C, VR_C, KT1_C, KT2_C
! Note: VBIAS_C is also in the common block in CURRENT*.FOR but not used by FERMIR/L
! REAL*8 VBIAS_C

  EXTERNAL FERMIR ! Declare FERMIR as an external function

! Setup common block variables for the right contact
  VR_C = 0.0D0       ! Chemical potential of right contact (e.g., in Rydberg)
  TEMP_R = 300.0D0   ! Temperature of right contact in Kelvin
  KB_RY = 8.6170D-5 / 13.60580D0 ! Boltzmann const in Rydberg/K
  KT2_C = KB_RY * TEMP_R

  ENERGY_VAL = 0.01D0 ! Energy at which to calculate f_R(E) (Rydberg)

  FERMI_DIST_R = FERMIR(ENERGY_VAL)
  PRINT *, "Fermi-Dirac distribution at E=", ENERGY_VAL, " for right contact is:", FERMI_DIST_R

  ENERGY_VAL = VR_C ! At the Fermi level
  FERMI_DIST_R = FERMIR(ENERGY_VAL)
  PRINT *, "Fermi-Dirac distribution at E=VR for right contact is:", FERMI_DIST_R

! Zero temperature case
  KT2_C = 0.0D0
  ENERGY_VAL = VR_C - 0.01D0
  FERMI_DIST_R = FERMIR(ENERGY_VAL)
  PRINT *, "Fermi-Dirac (T=0) at E < VR for right contact is:", FERMI_DIST_R
  ENERGY_VAL = VR_C + 0.01D0
  FERMI_DIST_R = FERMIR(ENERGY_VAL)
  PRINT *, "Fermi-Dirac (T=0) at E > VR for right contact is:", FERMI_DIST_R

END PROGRAM TEST_FERMIR
```

## Dependencies and Interactions

*   **COMMON /VL/VL/VR/VR/KT1/KT1/KT2/KT2/:** The function relies on `VR` (chemical potential of the right contact) and `KT2` (thermal energy `kB*TR` for the right contact) being correctly set in this common block by a calling routine or a setup module.
*   **Numerical Stability:** The implementation includes a branch for `KT2.EQ.0.0d0` (zero temperature), returning a step function (1 if `E < VR`, 0 if `E >= VR`). For `KT2 > 0`, it uses two different but mathematically equivalent forms for `E < VR` and `E >= VR` to potentially improve numerical stability when `(E-VR)/KT2` is large and positive or large and negative, preventing overflow in `dexp`.
    *   If `E < VR`: `1.0d0/(1.0d0+dexp((E-VR)/KT2))` (standard form)
    *   If `E >= VR`: `dexp(-(E-VR)/KT2)/(1.0d0+dexp(-(E-VR)/KT2))`, which is equivalent to `1.0d0/(dexp((E-VR)/KT2)+1.0d0)`.
*   This function is a fundamental component used by various subroutines (`CURRENT.FOR`, `CURRENT_EP.FOR`, `SCURRENT.FOR`) that calculate thermoelectric transport properties.
```
