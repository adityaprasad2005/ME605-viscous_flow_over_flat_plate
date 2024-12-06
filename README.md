# Computational Solution of Viscous Flow over a Flat Plate

This project simulates the viscous flow over a flat plate using the Finite Difference Method (FDM) for three different numerical schemes: Explicit Euler, Implicit Euler, and Crank-Nicolson. The problem statement considers a steady, laminar, and incompressible flow with a large Reynolds number, allowing for the boundary layer approximation.

**Files**

* CrankNicolson_trial.m: Script implementing the Crank-Nicolson scheme.
* Explicit_trial.m: Script implementing the Explicit Euler scheme.
* Implicit_trial.m: Script implementing the Implicit Euler scheme.

**Results**

The expected results are:

* Velocity (u) and y-velocity (v) fields visualized as contour plots.
* Normalized velocity (F') plotted as a line function of the similarity variable (η) and compared with the Blasius solution.
* Boundary layer thickness variation plotted against the x-coordinate and compared with the Blasius flat plate boundary layer solution.
