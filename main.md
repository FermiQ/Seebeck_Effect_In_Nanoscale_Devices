# main.for

## Overview

The `main.for` file contains the main program `Seebeck`, which serves as the primary driver for calculating thermoelectric properties, particularly the Seebeck coefficient, under different physical assumptions. The program prompts the user to select a computation mode: spinless, spin-polarized, or including electron-phonon interactions. Based on the selection, it reads relevant input data (transmission spectra, coupling parameters), sets up temperature and voltage parameters, and then calls the appropriate calculation subroutines (`CURRENT`, `SCURRENT`, `CURRENT_EP`, `thermcd`). The results, typically Seebeck coefficients or ZT values as a function of temperature or bias, are written to output files.

## Key Components

*   **PROGRAM Seebeck:**
    *   The main executable program.
    *   Initializes parameters, reads user input for calculation mode, reads data files, iterates over temperature or voltage, calls physics subroutines, and writes output.

## Important Variables/Constants

*   **L (Parameter):** Integer, set to 201. Defines the size of arrays holding energy grid and transmission data.
*   **NMODE (Parameter):** Integer, set to 3. Defines the number of phonon modes for electron-phonon interaction calculations.
*   **CASES (Input):** Integer. User input to select calculation mode:
    *   1: Spinless calculation (uses `CURRENT`).
    *   2: Spin-polarized calculation (uses `SCURRENT`).
    *   3: Electron-phonon interaction calculation, spinless (uses `CURRENT_EP` and `thermcd`).
*   **COMMON /VL/VL/VR/VR/KT1/KT1/KT2/KT2/VBIAS/VBIAS:**
    *   **VL, VR:** Real*8, chemical potentials of left/right contacts. Initialized based on `VBIAS`.
    *   **KT1, KT2:** Real*8, `kB*T` for left/right contacts.
    *   **VBIAS:** Real*8, applied bias voltage (in Volts, then converted or used in Rydberg context).
*   **x1(L):** Real*8 array, energy grid points (read from file).
*   **x4(L):** Real*8 array, transmission data for spinless case (read from `abscurrent.s1`).
*   **x4up(L), x4dn(L):** Real*8 arrays, transmission data for spin-up and spin-down channels (read from `46AP.txt`).
*   **CLR(NMODE), CRR(NMODE):** Real*8 arrays, electron-phonon coupling strengths for left/right (read from `couple.dat`).
*   **fr(NMODE):** Real*8 array, phonon frequencies (read from `w.input`, then converted to Rydberg).
*   **TL, TR:** Real*8, temperatures of left and right contacts (Kelvin). Varied in loops.
*   **Tw:** Real*8, phonon temperature (used in electron-phonon case, seems to be set to TR in the loop, but then `CURRENT_EP` is called with `Tw` as an argument which is not explicitly set in Case 3 loop before call, this might be a bug or `Tw` is intended to be passed from a prior calculation not shown). In `CASE(3)` loop, `Tw` is not explicitly assigned before calling `CURRENT_EP`. It will use whatever value `Tw` had from a previous context or if it was initialized (not shown here). This is a potential issue.
*   **SEE, SEE1, SEE2:** Real*8, output Seebeck coefficients from `CURRENT`, `SCURRENT`, or `CURRENT_EP`.
*   **KEQ, KQ, KQQ:** Real*8, thermal conductance related outputs.
*   **ZT:** Real*8, figure of merit `S^2 * G * T / K`, where G is electrical conductance (related to HH), K is thermal conductance.

## File I/O

*   **Input Files:**
    *   Reads calculation mode from user console.
    *   `abscurrent.s1`: Reads energy grid and transmission data (for spinless and EP cases). Columns: x1 (energy), x2, x3, x4 (transmission), x5, x6, x7.
    *   `46AP.txt`: Reads energy grid, spin-down transmission, spin-up transmission (for spin case).
    *   `couple.dat`: Reads electron-phonon coupling strengths `CLR`, `CRR` (for EP case).
    *   `w.input`: Reads phonon frequencies (for EP case, converted from cm^-1).
*   **Output Files:**
    *   `1.dat`: Writes energy (shifted by VL) and transmission (x4) - seems like a debug/check output.
    *   `SEE.dat`: Outputs `TR`, `SEE` for spinless case.
    *   `SEE_46AP.dat`: Outputs `TR`, `SEE`, `SEE1`, `SEE2` for spin case. Also prints to console.
    *   `ZT.dat`: Outputs `TR`, `ZT` for electron-phonon case.
    *   Prints "Task done!" or "TASK DONE" to console upon completion of loops.

## Usage Examples

The program is executed, and it interactively asks for the computation mode:
1.  User enters `1`, `2`, or `3`.
2.  The program reads necessary input files based on the choice.
3.  It then iterates through a temperature range (for modes 1 and 2) or a voltage range (for mode 3), calculating and writing thermoelectric properties to respective `.dat` files.

```
(Execution from a terminal)
> ./seebeck.out  (or equivalent executable name)
 BITTE WAEHLEN COMPUTATION MODE
 SPINLESS=(1), SPIN=(2), EP=(3)
1
 NOW SPINLESS CALCULATION
 (Program runs and creates SEE.dat)
 Task done!
```

## Dependencies and Interactions

*   **CURRENT (SUBROUTINE):** Called in `CASE(1)` for spinless calculations.
*   **SCURRENT (SUBROUTINE):** Called in `CASE(2)` for spin-polarized calculations.
*   **CURRENT_EP (SUBROUTINE):** Called in `CASE(3)` for calculations including electron-phonon interactions.
*   **thermcd (SUBROUTINE):** Called in `CASE(3)` to calculate `KQQ` (likely phonon thermal conductance).
*   **Common Block `/VL/VR/KT1/KT2/VBIAS/`:** Used to share contact properties (chemical potentials, temperatures, bias) with the calculation subroutines.
*   The program structure uses `SELECT CASE` for different physics models and `DO WHILE` loops to vary parameters like temperature or bias.
*   **Potential issue:** In `CASE(3)`, the variable `Tw` (phonon temperature) is passed to `CURRENT_EP` but does not appear to be explicitly set within that loop or case statement before the call. It relies on its value from outside the loop or a default initialization. The variable `NCH2` and `scaleM` passed to `thermcd` are also not initialized within `CASE(3)`.
```
