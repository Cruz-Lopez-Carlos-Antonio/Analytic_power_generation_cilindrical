# Analytical Solution for Power Generation in Cylindrical Microchannels with High Surface Zeta Potential

## Overview of the Repository

The present repository contains the **Python** implementation of the semianalytical solution developed for the electrohydrodynamic model of power generation in cylindrical microchannels with high surface zeta potential.

The computational procedure follows the analytical formulation step by step. The electric potential distribution within the electrical double layer is first obtained from the nonlinear Poisson–Boltzmann equation. Once this potential is known, the remaining hydrodynamic and electrokinetic quantities are reconstructed through a sequence of numerical quadratures, including the axial velocity profile, streaming potential, pressure distribution, and energy-conversion efficiency.

The implementation is primarily based on **NumPy** and **SciPy** as well as **Wolfram Mathematica** and is organized to allow the numerical verification of each stage of the semianalytical solution. In particular, the code includes consistency checks for the governing equations, boundary conditions, integral constraints, and numerical convergence.
