# FERMIL.FOR

## Overview

The `FERMIL.FOR` file defines the `FERMIL` function, which calculates the Fermi-Dirac distribution function for the left electron reservoir/contact. This function is analogous to `FERMIR.FOR` but uses the parameters for the left contact (`VL` and `KT1 = kB*TL`). The Fermi-Dirac distribution gives the probability that an electron state at a given energy `E` is occupied.

## Key Components

*   **REAL*8 FUNCTION FERMIL(E):**
    *   Calculates the Fermi-Dirac distribution `f_L(E) = 1 / (1 + exp((E - VL) / KT1))`.
    *   **E (Input):** Real*8, the energy at which to evaluate the distribution.
    *   **Returns (FERMIL):** Real*8, the value of the Fermi-Dirac distribution function for the left contact.

## Important Variables/Constants

*   **COMMON /VL/VL/VR/VR/KT1/KT1/KT2/KT2/:**
    *   **VL (Input, from common block):** Real*8, the chemical potential (Fermi level) of the left contact.
    *   **KT1 (Input, from common block):** Real*8, the thermal energy of the left contact (`kB * TL`, where `kB` is the Boltzmann constant and `TL` is the temperature of the left contact). This value is expected to be in the same energy units as `E` and `VL`.

## Usage Examples

The `FERMIL` function is typically called within integration loops in routines like `CURRENT.FOR`, `CURRENT_EP.FOR`, and `SCURRENT.FOR` to determine the occupation probability of states at different energies in the left contact.

```fortran
PROGRAM TEST_FERMIL
  IMPLICIT NONE
  REAL*8 ENERGY_VAL, FERMI_DIST_L

  COMMON /VL/VL_C/VR/VR_C/KT1_C/KT1_C/KT2_C/KT2_C
  REAL*8 VL_C, VR_C, KT1_C, KT2_C

  EXTERNAL FERMIL ! Declare FERMIL as an external function

! Setup common block variables for the left contact
  VL_C = 0.05D0    ! Chemical potential of left contact (e.g., in Rydberg)
  TEMP_L = 300.0D0 ! Temperature of left contact in Kelvin
  KB_RY = 8.6170D-5 / 13.60580D0 ! Boltzmann const in Rydberg/K
  KT1_C = KB_RY * TEMP_L

  ENERGY_VAL = 0.06D0 ! Energy at which to calculate f_L(E) (Rydberg)

  FERMI_DIST_L = FERMIL(ENERGY_VAL)
  PRINT *, "Fermi-Dirac distribution at E=", ENERGY_VAL, " for left contact is:", FERMI_DIST_L

  ENERGY_VAL = VL_C ! At the Fermi level
  FERMI_DIST_L = FERMIL(ENERGY_VAL)
  PRINT *, "Fermi-Dirac distribution at E=VL for left contact is:", FERMI_DIST_L

! Zero temperature case
  KT1_C = 0.0D0
  ENERGY_VAL = VL_C - 0.01D0
  FERMI_DIST_L = FERMIL(ENERGY_VAL)
  PRINT *, "Fermi-Dirac (T=0) at E < VL for left contact is:", FERMI_DIST_L
  ENERGY_VAL = VL_C + 0.01D0
  FERMI_DIST_L = FERMIL(ENERGY_VAL)
  PRINT *, "Fermi-Dirac (T=0) at E > VL for left contact is:", FERMI_DIST_L

END PROGRAM TEST_FERMIL
```

## Dependencies and Interactions

*   **COMMON /VL/VL/VR/VR/KT1/KT1/KT2/KT2/:** The function relies on `VL` (chemical potential of the left contact) and `KT1` (thermal energy `kB*TL` for the left contact) being correctly set in this common block.
*   **Numerical Stability:** Similar to `FERMIR.FOR`, this function implements logic to handle the `KT1.EQ.0.0d0` (zero temperature) case, returning a step function. For `KT1 > 0`, it uses two different but mathematically equivalent forms for `E < VL` and `E >= VL` to potentially improve numerical stability and prevent overflow in `dexp`.
    *   If `E < VL`: `1.0d0/(1.0d0+dexp((E-VL)/KT1))`
    *   If `E >= VL`: `dexp(-(E-VL)/KT1)/(1.0d0+dexp(-(E-VL)/KT1))`
*   This function is a fundamental component used by various subroutines (`CURRENT.FOR`, `CURRENT_EP.FOR`, `SCURRENT.FOR`) that calculate thermoelectric transport properties.
```
