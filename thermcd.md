# thermcd.for

## Overview

The `thermcd.for` file defines the `thermcd` subroutine and the `FCR` function, as well as an integer function `MACHPD`. The `thermcd` subroutine appears to calculate thermal conductance or related properties, possibly involving phonon contributions and surface effects, given parameters like temperature of left/right contacts (TL, TR) and material properties. It makes use of the `ZEROIN` subroutine for root finding and `SURFPH` for surface phonon calculations.

## Key Components

*   **SUBROUTINE thermcd(NCH2, TL, TR, scaleM, KQQ):**
    *   Calculates thermal properties.
    *   **NCH2:** An integer input, possibly related to the number of CH2 units or a geometric parameter of the nanoscale device.
    *   **TL, TR:** Double precision inputs representing temperatures of the left and right contacts/reservoirs.
    *   **scaleM:** A double precision input, likely a scaling factor for current or heating.
    *   **KQQ:** A double precision output, possibly the calculated thermal conductance or a related quantity.
*   **FUNCTION FCR(X):**
    *   A double precision function that defines a polynomial equation: `X**6 - 8*X**4 + 8*(3 - 2*V*V)*X*X - 16*(1 - V*V)`.
    *   This function is used as input to the `ZEROIN` subroutine to find its roots.
    *   **X:** Double precision input variable.
    *   Depends on `VL` and `VT` (longitudinal and transverse sound velocities) through the common block `/SURF/`.
*   **INTEGER FUNCTION MACHPD(X):**
    *   Determines if `1.0D0 < X`. Returns `1` if true, `0` otherwise. This is a simple machine precision check, likely used by `ZEROIN`.
    *   **X:** Double precision input.

## Important Variables/Constants

*   **COMMON/SURF/VL, VT, CR:**
    *   **VL:** Longitudinal sound velocity.
    *   **VT:** Transverse sound velocity.
    *   **CR:** Rayleigh wave velocity (calculated using `ZEROIN` on `FCR`).
*   **COMMON/TVSPOWER/vTw, vPTOE:**
    *   **vTw(2000):** Real*8 array, likely storing temperature values.
    *   **vPTOE(2000):** Real*8 array, likely storing power dissipated to electrode.
*   **AKB:** Boltzmann constant in units of (eV/K) / 13.6 (Hartree).
*   **hbar:** Planck's constant (set to 1.0d0, implying atomic units or normalization).
*   **PI:** Mathematical constant Pi.
*   **RHO:** Density of the material (e.g., 19.3d0 for Au, 2.7d0 for Al).
*   **Ymdl:** Young's modulus or a related elastic constant.
*   **dcross, dlong:** Dimensions, likely cross-sectional area and length of the device.
*   **Cth:** A thermal capacitance or related coefficient calculated using G1, G2, G3 (from `SURFPH`), CR, VL, VT, RHO, and DHBAR.
*   **KB:** Boltzmann constant in CGS units (1.38065D-16 erg/K).
*   **DHBAR:** Planck's constant in CGS units (1.05459d-27 erg*s).
*   **KK:** A constant, possibly related to spring constant or coupling, calculated using `dcross`, `dlong`, and `Ymdl`.
*   **betaM:** A scaling factor (0.78d0 * 0.529d0), potentially related to contact geometry or decay length from Mark Reed PRB 68 035416 (2003).
*   **scaleM:** Input scaling factor, described in comments as `DEXP(-1.d0*betaM*(ddd2-ddd1))`, used to scale current and heating.

## Usage Examples

The `thermcd` subroutine would typically be called from a main program that sets up the physical parameters of the system.

```fortran
PROGRAM MAIN_CALC
  IMPLICIT DOUBLE PRECISION (A-H,O-Z)
  IMPLICIT INTEGER (I-N)
  REAL*8 TEMPERATURE_L, TEMPERATURE_R, SCALING_M, THERMAL_COND_Q
  INTEGER N_CH2_UNITS

  N_CH2_UNITS = 4
  TEMPERATURE_L = 300.0D0 ! Kelvin
  TEMPERATURE_R = 310.0D0 ! Kelvin
  SCALING_M = 1.0D0     ! No additional scaling

  CALL thermcd(N_CH2_UNITS, TEMPERATURE_L, TEMPERATURE_R, SCALING_M, THERMAL_COND_Q)

  PRINT *, "Calculated KQQ:", THERMAL_COND_Q

  STOP
END PROGRAM MAIN_CALC
```

## Dependencies and Interactions

*   **ZEROIN (SUBROUTINE):** Called to find the root of the `FCR` function, which determines the Rayleigh wave velocity `CR`.
*   **SURFPH (SUBROUTINE):** Called to calculate G1, G2, G3, which are likely related to surface phonon density of states or dispersion. (Content of `SURFPH.F90` would clarify).
*   **F_R, F_L, FCR, FCNPOWER2, FCNPOWERnod (EXTERNAL REAL*8 FUNCTIONS):** Declared as external, suggesting they are provided elsewhere and used by `thermcd` or its callees (though not directly visible in the provided `thermcd` code snippet, they might be passed to other routines or used in parts of the code not shown). However, `FCR` is defined within this file. The other F_ functions are not directly used in `thermcd` but might be part of the larger context this common block is used in.
*   **Common Blocks:**
    *   `/SURF/`: Shares `VL`, `VT`, `CR` with other routines (e.g., `FCR` function within this file, and likely `SURFPH`).
    *   `/TVSPOWER/`: Shares `vTw`, `vPTOE` arrays, suggesting these are populated or used by other parts of the simulation to track temperature and power dissipation.
*   The subroutine sets material parameters (sound velocities, density) for Au or Al (commented out), indicating it's configurable for different materials.
```
