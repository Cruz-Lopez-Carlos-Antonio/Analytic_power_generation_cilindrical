# Analytical Solution for Power Generation in Cylindrical Microchannels with High Surface Zeta Potential

## Overview of the Repository

The present repository contains the **Python** and **Mathematica** implementations of the semianalytical solution developed for the electrohydrodynamic model of power generation in cylindrical microchannels with high surface zeta potential.

The computational procedure follows the analytical formulation step by step. The electric potential distribution within the electrical double layer is first obtained from the nonlinear Poisson–Boltzmann equation. Once this potential is known, the remaining hydrodynamic and electrokinetic quantities are reconstructed through a sequence of numerical quadratures, including the axial velocity profile, streaming potential, pressure distribution, and energy-conversion efficiency.

The implementation is primarily based on **NumPy** and **SciPy** as well as **Wolfram Mathematica** and is organized to allow the numerical verification of each stage of the semianalytical solution. In particular, the code includes consistency checks for the governing equations, boundary conditions, integral constraints, and numerical convergence.
## Theoretical Framework

The computational implementation is based on the dimensionless formulation of the electrohydrodynamic transport problem in a cylindrical microchannel. The model couples the electric potential distribution within the electrical double layer with the hydrodynamic and electrokinetic fields responsible for pressure-driven flow and streaming-potential generation.

### Physical Parameters

The physical parameters employed in the numerical implementation are summarized below.

| Parameter | Symbol | Value | Unit |
|---|---:|---:|---|
| Microchannel length | $l$ | $1.0\times10^{-3}$ | $\mathrm{m}$ |
| Microchannel radius | $a$ | $1.0\times10^{-8}$ | $\mathrm{m}$ |
| Temperature | $T$ | $298$ | $\mathrm{K}$ |
| Fluid density | $\rho$ | $1000$ | $\mathrm{kg\,m^{-3}}$ |
| Dynamic viscosity | $\mu$ | $0.891\times10^{-3}$ | $\mathrm{Pa\,s}$ |
| Diffusion coefficient | $D$ | $1.312\times10^{-9}$ | $\mathrm{m^2\,s^{-1}}$ |
| Permittivity | $\varepsilon$ | $6.954\times10^{-10}$ | $\mathrm{C\,V^{-1}\,m^{-1}}$ |
| Ionic valence | $z$ | $1$ | — |
| Bulk osmotic concentration | $C_{\mathrm{osD},b}$ | $0.5$ | $\mathrm{M}$ |
| Water flux | $J_w$ | $8.68327\times10^{-7}$ | $\mathrm{m\,s^{-1}}$ |
| pH | $\mathrm{pH}$ | $7.0$ | — |
| Dissociation constant | $\mathrm{p}K$ | $7.5$ | — |
| Surface site density | $\Gamma$ | $8$ | $\mathrm{nm^{-2}}$ |
| Stern-layer capacitance | $C_{\mathrm{Stern}}$ | $2.9$ | $\mathrm{F\,m^{-2}}$ |
| Zeta potential | $\zeta$ | $-0.03165$ | $\mathrm{V}$ |
| Viscoelectric coefficient | $f$ | $2.3\times10^{-16}$ | $\mathrm{m^2\,V^{-2}}$ |

These physical quantities are used directly to reconstruct the dimensionless parameters required by the governing equations. In particular, the dimensionless quantities are **not introduced independently as numerical input values**, but are computed from the corresponding dimensional properties.

### Dimensionless Parameters

For the physical conditions considered in the implementation, the reference dimensionless parameters are

| Parameter | Symbol | Reference value |
|---|---:|---:|
| Reynolds number | $\mathrm{Re}$ | $9.7455\times10^{-4}$ |
| Aspect ratio | $\xi$ | $1.0\times10^{-5}$ |
| Dimensionless Debye parameter | $\delta$ | $23.2419$ |
| Pressure-related parameter | $\alpha$ | $3.2025\times10^{12}$ |
| Electrohydrodynamic coupling parameter | $\Lambda$ | $0.759514$ |
| Viscoelectric parameter | $\omega$ | $0.00151706$ |

In particular, the geometrical aspect ratio is

$$
\xi=\frac{a}{l},
$$

whereas

$$
\delta=a\kappa,
$$

with $\kappa^{-1}$ denoting the Debye length. Thus, $\delta$ represents the ratio between the microchannel radius and the characteristic Debye length.

The dimensionless surface potential is defined as

$$
\Psi_s=\frac{\zeta}{\zeta_T},
$$

where the thermal potential is

$$
\zeta_T=\frac{k_B T}{ze}.
$$

Here, $k_B$ is the Boltzmann constant and $e$ is the elementary charge. The notation $\Psi_s$ is used for the dimensionless surface potential to distinguish it from the dimensionless axial coordinate $Z$.

### Dimensionless Governing Equations

The electric potential within the electrical double layer is governed by the cylindrical Poisson–Boltzmann equation

$$
\frac{d^2\Psi}{dR^2}
+
\frac{1}{R}\frac{d\Psi}{dR}
=
\delta^2\sinh\!\left(\Psi\right),
\qquad 0\leq R\leq1,
$$

subject to the boundary conditions

$$
\left.\frac{d\Psi}{dR}\right|_{R=0}=0,
\qquad
\Psi(1)=\Psi_s.
$$

The dimensionless viscoelectric correction to the viscosity is represented by

$$
M(R)
=
\exp\!\left[
\omega
\left(
\frac{d\Psi}{dR}
\right)^2
\right].
$$

Once the electrostatic potential is known, the function $F(R)$ governing the axial velocity field satisfies

$$
0=
-\Pi_D
+
\frac{1}{R}
\frac{d}{dR}
\left[
M(R)R\frac{dF}{dR}
\right]
-
\Lambda\delta^2
\sinh\!\left(\Psi(R)\right)\Omega,
$$

where

$$
\Pi_D=\alpha\xi^2,
$$

and the electrohydrodynamic coupling quantity $\Omega$ is defined through the integral relation

$$
\Omega
=
\int_0^1
F(R)\sinh\!\left(\Psi(R)\right)R\,dR.
$$

The corresponding boundary conditions for $F(R)$ are

$$
\left.\frac{dF}{dR}\right|_{R=0}=0,
\qquad
F(1)=0.
$$

#### Poisson–Boltzmann Equation

The electric potential within the electrical double layer is governed by the cylindrical Poisson–Boltzmann equation

$$\frac{d^2\Psi}{dR^2}+\frac{1}{R}\frac{d\Psi}{dR}=\delta^2\sinh\left(\Psi\right), \qquad 0\leq R\leq1.$$

subject to the boundary conditions

$$\left.\frac{d\Psi}{dR}\right|_{R=0}=0, \qquad \Psi(1)=\Psi_s.$$

#### Viscoelectric Correction

The dimensionless viscoelectric correction to the viscosity is represented by

$$M(R)=\exp\left[\omega\left(\frac{d\Psi}{dR}\right)^2\right].$$

#### Axial Flow Equation

Once the electrostatic potential is known, the function $F(R)$ governing the axial velocity field satisfies

$$0=-\Pi_D+\frac{1}{R}\frac{d}{dR}\left[M(R)R\frac{dF}{dR}\right]-\Lambda\delta^2\sinh\left(\Psi(R)\right)\Omega,$$

where

$$\Pi_D=\alpha\xi^2,$$

and the electrohydrodynamic coupling quantity $\Omega$ is defined through the integral relation

$$\Omega=\int_0^1F(R)\sinh\left(\Psi(R)\right)R\,dR.$$

The corresponding boundary conditions for $F(R)$ are

$$\left.\frac{dF}{dR}\right|_{R=0}=0, \qquad F(1)=0.$$

## Semianalytical Solution and Computational Implementation

Once the dimensionless electric potential $\Psi(R)$ has been obtained from the Poisson–Boltzmann equation, the remaining electrohydrodynamic problem can be reconstructed through a sequence of numerical quadratures.

First, the auxiliary function

$$H(R)=\int_0^R \tau\sinh\left(\Psi(\tau)\right)\,d\tau$$

is introduced. The solution for $F(R)$ can then be decomposed as

$$F(R)=F_0(R)+\Omega F_1(R),$$

where

$$F_0(R)=-\int_R^1\frac{\Pi_D s}{2M(s)}\,ds,$$

and

$$F_1(R)=-\Lambda\delta^2\int_R^1\frac{H(s)}{sM(s)}\,ds.$$

To determine the global coupling quantity $\Omega$, the auxiliary integrals

$$A=\int_0^1F_0(R)\sinh\left(\Psi(R)\right)R\,dR,$$

and

$$B=\int_0^1F_1(R)\sinh\left(\Psi(R)\right)R\,dR$$

are evaluated. Substitution into the integral definition of $\Omega$ leads to the compact relation

$$\Omega=\frac{A}{1-B},$$

provided that $1-B\neq0$. The function $F(R)$ is subsequently reconstructed from $F_0(R)$, $F_1(R)$, and $\Omega$.

### Pressure Gradient and Velocity Field

The dimensionless modified-pressure gradient is determined from the flow-rate constraint

$$\frac{1}{2}=\frac{d\widetilde{\Pi}}{dZ}\int_0^1RF(R)\,dR.$$

Defining

$$G:=\frac{d\widetilde{\Pi}}{dZ},$$

the pressure gradient can therefore be written as

$$G=\frac{1}{2\displaystyle\int_0^1RF(R)\,dR}.$$

The axial velocity profile is then reconstructed as

$$V_Z(R)=GF(R),$$

which satisfies the dimensionless flow-rate condition

$$\int_0^1RV_Z(R)\,dR=\frac{1}{2}.$$

Since $G$ is constant, the axial velocity is independent of the axial coordinate $Z$. The continuity equation, together with the impermeability condition at the microchannel wall, consequently gives

$$V_R(R,Z)=0.$$

### Streaming Potential and Pressure Distribution

The dimensionless streaming potential satisfies

$$\frac{d\Phi}{dZ}=-\Lambda\delta^2G\Omega.$$

Using the reference condition $\Phi(0)=0$, its axial distribution becomes

$$\Phi(Z)=-\Lambda\delta^2G\Omega Z.$$

Similarly, imposing $\widetilde{\Pi}(0)=1$ gives the modified pressure distribution

$$\widetilde{\Pi}(Z)=GZ+1.$$

The total dimensionless pressure is finally recovered as

$$\Pi(R,Z)=\widetilde{\Pi}(Z)+\cosh\left(\Psi(R)\right).$$

### Energy-Conversion Efficiency

The dimensionless energy-conversion efficiency is evaluated from

$$\eta=\eta_c\frac{\left(d\Phi/dZ\right)^2}{-d\Pi/dZ}.$$

Since $\cosh\left(\Psi(R)\right)$ is independent of $Z$, it follows that $d\Pi/dZ=G$, and therefore

$$\eta=\eta_c\frac{\left(\Lambda\delta^2G\Omega\right)^2}{-G}.$$

The parameter $\eta_c$ is reconstructed from the physical properties employed in the model rather than introduced independently as a numerical input.
