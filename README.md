# Analytical Solution for Power Generation in Cylindrical Microchannels with High Surface Zeta Potential

## Overview of the Repository

The present repository contains the **Python** and **Mathematica** implementations of the semianalytical solution developed for the electrohydrodynamic model of power generation in cylindrical microchannels with high surface zeta potential.

The computational procedure follows the analytical formulation step by step. The electric potential distribution within the electrical double layer is first obtained from the nonlinear Poisson–Boltzmann equation. Once this potential is known, the remaining hydrodynamic and electrokinetic quantities are reconstructed through a sequence of numerical quadratures, including the axial velocity profile, streaming potential, pressure distribution, and energy-conversion efficiency.

The implementation is primarily based on **NumPy** and **SciPy** as well as **Wolfram Mathematica** and is organized to allow the numerical verification of each stage of the semianalytical solution. In particular, the code includes consistency checks for the governing equations, boundary conditions, integral constraints, and numerical convergence.
## Theoretical Framework

The computational implementation is based on the dimensionless formulation of the electrohydrodynamic transport problem in a cylindrical microchannel. The model couples the electric potential distribution within the electrical double layer with the hydrodynamic and electrokinetic fields responsible for pressure-driven flow and streaming-potential generation.

## Overview of the Repository

<div style="background:#f1f7ff; padding:12px 15px; border-left:4px solid #4a90e2; border-radius:8px; margin-top:20px;">
  <strong>📘 Full Documentation</strong><br>
  The complete project website (opens in a new page) is available here:<br><br>
  👉 <a href="https://cruz-lopez-carlos-antonio.github.io/Analytic_power_generation_cilindrical/"
        target="_blank" rel="noopener noreferrer">
        https://cruz-lopez-carlos-antonio.github.io/Analytic_power_generation_cilindrical/
      </a>
</div>
