---
layout: default
title: Semi-analytic Solution
math: true
---

# Semi-analytic Solution

This section outlines the semianalytical approach used to calculate the parametric average streaming potential ($\Phi_{\mathrm{Av}}$) and the conversion efficiency ($\eta$).

<div style="background:#fff3cd; padding:15px; border-left:4px solid #ffeeba; border-radius:6px; margin:20px 0; color:#856404;">
  <strong>⚠️ Note:</strong><br>
  The following equations summarize the core of the semi-analytic reduction. For the full theoretical derivation, auxiliary demonstrations, and a detailed discussion on the separation of variables, please consult the original manuscript.
</div>

## The Fully Developed Assumption

The mathematical reduction begins by assuming a hydrodynamically fully developed regime. Under this assumption, the axial velocity profile remains constant along the longitudinal direction ($\partial V_Z/\partial Z=0$), which directly implies that the radial velocity component vanishes completely ($V_R=0$).

## Reduced Equation System

Applying the fully developed conditions, the original system can be reduced to a single Fredholm ordinary integro-differential equation for the function $F(R)$:

<div style="background:#f7f7f7; padding:15px; border-left:4px solid #4a90e2; border-radius:6px; margin:20px 0;">
$$
0=-\Pi_D +\frac{1}{R}\frac{d}{dR} \left(\exp\!\left[\omega \left(\Psi'(R)\right)^2\right] R\frac{dF}{dR} \right)- \Lambda\delta^2\sinh(\Psi) \int_0^1F(s) \sinh(\Psi(s))s\,ds
$$
</div>

This integro-differential equation is simplified into a pure-differential one by introducing the constant parameter $\Omega$, defined as:[cite: 7]

<div style="background:#f7f7f7; padding:15px; border-left:4px solid #4a90e2; border-radius:6px; margin:20px 0;">
$$
\Omega=\int_0^1F(R) \sinh(\Psi(R))R\,dR
$$
</div>

## Auxiliary Functions

To facilitate the analytical integration, two auxiliary functions are defined for the integral of the electrical forces and the viscoelectric modification, respectively:[cite: 7]

<div style="background:#f7f7f7; padding:15px; border-left:4px solid #4a90e2; border-radius:6px; margin:20px 0;">
$$
H(R)=\int_{0}^{R}\tau \sinh(\Psi(\tau))d\tau, \qquad M(R)=\exp(\omega[\Psi'(R)]^2)
$$
</div>

## Solution for $F(R)$

Integrating the modified momentum equation and applying the no-slip boundary condition at the wall ($F(1)=0$), the radial function $F(R)$ is obtained as the superposition of two components:[cite: 7]

<div style="background:#f7f7f7; padding:15px; border-left:4px solid #4a90e2; border-radius:6px; margin:20px 0;">
$$
F(R)=F_0(R)+\Omega F_1(R)
$$
</div>

where $F_0(R)$ and $F_1(R)$ are given by:[cite: 7]

<div style="background:#f7f7f7; padding:15px; border-left:4px solid #4a90e2; border-radius:6px; margin:20px 0;">
$$
F_0(R)=-\int_{R}^{1}\frac{\Pi_D x}{2M(x)}dx \quad\text{and} \quad F_1(R)=-\Lambda\delta^2\int_{R}^{1}\frac{H(x)}{xM(x)}dx
$$
</div>

## Streaming Potential and Conversion Efficiency

The average streaming potential can be evaluated directly using the previous definitions:[cite: 7]

<div style="background:#f7f7f7; padding:15px; border-left:4px solid #4a90e2; border-radius:6px; margin:20px 0;">
$$
\Phi_{\mathrm{Av}} =-\frac{\Lambda\delta^2}{2} \frac{d\widetilde{\Pi}}{dZ}\frac{A}{1-B}
$$
</div>

where the constants $A$ and $B$ encapsulate the integrals of the $F_0$ and $F_1$ functions across the microchannel section:[cite: 7]

<div style="background:#f7f7f7; padding:15px; border-left:4px solid #4a90e2; border-radius:6px; margin:20px 0;">
$$
A=\int_{0}^{1}F_0(R)\sinh\left(\Psi\left(R\right)\right)R\,dR \quad\text{and}\quad B=\int_{0}^{1}F_1(R)\sinh\left(\Psi\left(R\right)\right)R\,dR
$$
</div>

Finally, the conversion efficiency $\eta$, representing the ratio of output electrical power to input hydrodynamic power, is formulated as:[cite: 7]

<div style="background:#f7f7f7; padding:15px; border-left:4px solid #4a90e2; border-radius:6px; margin:20px 0;">
$$
\eta=\eta_c \frac{\left(\dfrac{d\Phi}{dZ} \right)^2}{\left( -\dfrac{d\widetilde{\Pi}}{dZ}\right)}
$$
</div>
