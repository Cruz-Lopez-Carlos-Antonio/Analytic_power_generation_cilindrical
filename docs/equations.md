---
layout: default
title: Physical Model
math: true
---
<div style="background:#fff3cd; padding:15px; border-left:4px solid #ffeeba; border-radius:6px; margin:20px 0; color:#856404;">
  <strong>⚠️ Disclaimer:</strong><br>
  The present repository focuses exclusively on the computational implementation and numerical validations of this project. The following sections list the primary mathematical identities governing the physical system. For a comprehensive discussion of the underlying theory, detailed derivations, and physical implications, the interested reader is strongly encouraged to consult the original manuscript, which was recently submitted to the journal <em>Physica Scripta</em>.
</div>

## Physical Model and Governing Equations

The physical system consists of a cylindrical microchannel of radius $a$ and length $l$. The inner surface of the microchannel is assumed to be negatively charged. As the electrolyte flows through the channel, the mobile charge within the electrical double layer is transported in the axial direction, leading to charge accumulation at the channel ends and thereby generating a streaming potential.

<div style="text-align: center; margin: 30px 0;">
  <!-- Asegúrate de subir tu diagrama compilado como imagen a la carpeta de assets -->
  <img src="{{ '/assets/images/diagrama_microcanal.png' | relative_url }}" alt="Schematic representation of the cylindrical microchannel" style="max-width: 100%; border-radius: 8px; box-shadow: 0 4px 6px rgba(0,0,0,0.05);">
  <p style="color: #666; font-size: 0.9rem; margin-top: 10px;"><em>Figure 1: Schematic representation of the cylindrical microchannel.</em></p>
</div>

## Original Mass Balance System

Before introducing the dimensionless variables, the fundamental physics of the cylindrical microchannel is governed by the original dimensional system of equations.[cite: 4] For a steady, incompressible, and axisymmetric flow, the continuity equation is given by:[cite: 4]

<div style="background:#f7f7f7; padding:15px; border-left:4px solid #4a90e2; border-radius:6px; margin:20px 0;">
$$
\frac{1}{r}\frac{\partial}{\partial r}\left(r v_r\right)+\frac{\partial v_z}{\partial z} = 0
$$
</div>

The electric double-layer (EDL) potential, which varies in the radial direction normal to the microchannel surface, is governed by the fully nonlinear Poisson–Boltzmann equation in cylindrical coordinates:[cite: 4]

<div style="background:#f7f7f7; padding:15px; border-left:4px solid #4a90e2; border-radius:6px; margin:20px 0;">
$$
\frac{d^2\psi}{dr^2} +\frac{1}{r}\frac{d\psi}{dr} =\frac{2zen_{\infty}}{\varepsilon} \sinh\!\left(\frac{ze\psi}{k_B T}\right)
$$
</div>

By incorporating the radial electrostatic contribution into a modified pressure $\widetilde{p}$, and substituting the constitutive relations for the viscous stresses and electric body-force components, the radial and axial momentum conservation equations are respectively formulated as:[cite: 4]

<div style="background:#f7f7f7; padding:15px; border-left:4px solid #4a90e2; border-radius:6px; margin:20px 0;">
$$
\begin{aligned}
\rho\left(v_r\frac{\partial v_r}{\partial r} +v_z\frac{\partial v_r}{\partial z}\right) &= -\frac{\partial\widetilde{p}}{\partial r} +\frac{2}{r}\frac{\partial}{\partial r} \left(\mu r\frac{\partial v_r}{\partial r} \right) \\
&\quad +\frac{\partial}{\partial z} \left[\mu\left(\frac{\partial v_r}{\partial z} +\frac{\partial v_z}{\partial r}\right) \right]-2\mu\frac{v_r}{r^2}
\end{aligned}
$$
</div>

<div style="background:#f7f7f7; padding:15px; border-left:4px solid #4a90e2; border-radius:6px; margin:20px 0;">
$$
\begin{aligned}
\rho\left( v_r\frac{\partial v_z}{\partial r} +v_z\frac{\partial v_z}{\partial z} \right) &= -\frac{\partial\widetilde{p}}{\partial z} +\frac{1}{r}\frac{\partial}{\partial r} \left[\mu r\left( \frac{\partial v_r}{\partial z} +\frac{\partial v_z}{\partial r} \right)\right] \\
&\quad + 2\frac{\partial}{\partial z} \left(\mu\frac{\partial v_z}{\partial z}\right)+ 2zen_{\infty} \sinh\!\left( \frac{ze\psi}{k_B T} \right) \frac{\partial\phi}{\partial z}
\end{aligned}
$$
</div>


### Dimensionless Variables

To generalize the solution, the governing equations are transformed into a dimensionless form. The mathematical formulation uses the dimensionless variables defined in the following table:

| Dimensionless variable | Definition | Description |
| :--- | :--- | :--- |
| $Z$ | $z/l$ | Where $z$ is the axial coordinate and $l$ is the characteristic channel length. |
| $R$ | $r/a$ | Where $r$ is the radial coordinate and $a$ is the radius of the microchannel.|
| $V_Z$ | $v_z/J_w$ | Where $v_z$ is the axial velocity component and $J_w$ is defined via the volumetric flow rate as $Q=A_m J_w$.|
| $V_R$ | $v_r/v_{r,c}$ | Where $v_r$ is the radial velocity component and $v_{r,c}$ is the characteristic radial velocity scale.|
| $\Psi$ | $\psi/\zeta_T$ | Where $\psi$ is the electric double-layer (EDL) potential and $\zeta_T$ is the thermal potential.|
| $\Phi$ | $\phi/\phi_c$ | Where $\phi$ is the streaming potential and $\phi_c$ is the characteristic streaming potential scale.|
| $\widetilde{\Pi}$ | $\widetilde{p}/P_{\mathrm{OS,D,b}}$ | Where $\widetilde{p}$ is the modified pressure and $P_{\mathrm{OS,D,b}}$ is the bulk osmotic pressure reference scale.|

### Dimensionless Balance Equations

The solution for the streaming potential, the velocity field, and the induced pressure is obtained by solving the following system of differential equations:

<div style="background:#f7f7f7; padding:15px; border-left:4px solid #4a90e2; border-radius:6px; margin:20px 0;">
$$
\frac{d^2\Psi}{dR^2} + \frac{1}{R}\frac{d\Psi}{dR} = \delta^2\sinh(\Psi)
$$
</div>

<div style="background:#f7f7f7; padding:15px; border-left:4px solid #4a90e2; border-radius:6px; margin:20px 0;">
$$
\frac{1}{R}\frac{\partial}{\partial R}(RV_R) + \frac{\partial V_Z}{\partial Z} = 0
$$
</div>

<div style="background:#f7f7f7; padding:15px; border-left:4px solid #4a90e2; border-radius:6px; margin:20px 0;">
$$
\begin{aligned}
\mathrm{Re}\,\xi^2 \left(V_R\frac{\partial V_R}{\partial R} +V_Z\frac{\partial V_R}{\partial Z}\right) 
&= -\alpha \frac{\partial\widetilde{\Pi}}{\partial R} +\frac{2}{R} \frac{\partial}{\partial R}\left( \exp\!\left[\omega \left(\frac{d\Psi}{dR}\right)^2 \right]R\frac{\partial V_R}{\partial R} \right) \\
&\quad + \frac{\partial}{\partial Z} \left(\exp\!\left[\omega\left( \frac{d\Psi}{dR}\right)^2\right] \left[\xi^2\frac{\partial V_R}{\partial Z} +\frac{\partial V_Z}{\partial R}\right] \right) \\
&\quad - 2\exp\!\left[\omega\left(\frac{d\Psi}{dR} \right)^2\right] \frac{V_R}{R^2}
\end{aligned}
$$
</div>

<div style="background:#f7f7f7; padding:15px; border-left:4px solid #4a90e2; border-radius:6px; margin:20px 0;">
$$
\begin{aligned}
\mathrm{Re}\,\xi^2 \left( V_R\frac{\partial V_Z}{\partial R} + V_Z\frac{\partial V_Z}{\partial Z} \right)
&= -\alpha\xi^2 \frac{\partial\widetilde{\Pi}}{\partial Z} + \frac{1}{R} \frac{\partial}{\partial R} \left( \exp\!\left[ \omega \left( \frac{d\Psi}{dR} \right)^2 \right] R \left[ \xi^2 \frac{\partial V_R}{\partial Z} + \frac{\partial V_Z}{\partial R} \right] \right) \\
&\quad + 2\xi^2 \frac{\partial}{\partial Z} \left( \exp\!\left[ \omega \left( \frac{d\Psi}{dR} \right)^2 \right] \frac{\partial V_Z}{\partial Z} \right) + \sinh(\Psi)\frac{d\Phi}{dZ}
\end{aligned}
$$
</div>

<div style="background:#f7f7f7; padding:15px; border-left:4px solid #4a90e2; border-radius:6px; margin:20px 0;">
$$
\frac{d\Phi}{dZ}=-\Lambda\delta^2 \int_0^1V_Z(R)\sinh(\Psi(R))R\,dR
$$
</div>

<div style="background:#f7f7f7; padding:15px; border-left:4px solid #4a90e2; border-radius:6px; margin:20px 0;">
$$
\frac{1}{2}=\int_{0}^{1}RV_Z(R)dR
$$
</div>

### Boundary Conditions

The set of governing differential equations is solved subject to the following boundary conditions:

<div style="background:#f7f7f7; padding:15px; border-left:4px solid #4a90e2; border-radius:6px; margin:20px 0;">
$$
\begin{aligned}
\Psi(1) &= \Psi_s &&\quad \text{and} \quad& \Psi'(0) &= 0, \\
V_Z(1,Z) &= 0 &&\quad \text{and} \quad& \partial_RV_Z(0,Z) &= 0, \\
V_R(1,Z) &= 0 &&\quad \text{and} \quad& \partial_RV_R(0,Z) &= 0, \\
\widetilde{\Pi}(R,0) &= 1 &&\quad \text{and} \quad& \Phi(0) &= 0.
\end{aligned}
$$
</div>

At the microchannel wall ($R=1$), the dimensionless EDL potential $\Psi$ is determined by the surface potential $\Psi_s$. Additionally, the axial and radial velocities ($V_Z$ and $V_R$) must satisfy the standard no-slip and no-penetration constraints at this boundary. Due to the axisymmetric geometry of the system, the radial derivatives of $\Psi$, $V_Z$, and $V_R$ all vanish at the central axis ($R=0$). Lastly, at the inlet of the electrokinetic region ($Z=0$), the modified pressure $\widetilde{\Pi}$ equals the bulk osmotic pressure, and the streaming potential $\Phi$ is set to zero, which reflects a state of zero net charge accumulation.

### Dimensionless Parameters

The mathematical formulation relies on several derived dimensionless parameters, whose definitions and numerical magnitudes for the current model are summarized in the table below:

| Parameter | Symbol | Definition | Value |
| :--- | :--- | :--- | :--- |
| Reynolds number | $\mathrm{Re}$ | $\frac{\rho J_w l}{\mu_0}$ | $9.7455 \times 10^{-4}$ |
| Aspect ratio | $\xi$ | $\frac{a}{l}$ | $1.0000 \times 10^{-5}$ |
| Electrokinetic radius | $\delta$ | $a\kappa$ | $2.3242 \times 10^{1}$ |
| Dimensionless reference pressure | $\alpha$ | $\frac{P_{\mathrm{OS,D,b}}\,l}{\mu_0 J_w}$ | $3.2025 \times 10^{12}$ |
| Electrokinetic coupling parameter | $\Lambda$ | $\frac{2\varepsilon^{2}\kappa^{2}\zeta_{T}^{2}}{\mu_{0}\sigma_{\infty}F_{cc}}$ | $7.5951 \times 10^{-1}$ |
| Viscoelectric parameter | $\omega$ | $f_{\mathrm{VE}} \left(\frac{\zeta_T}{a}\right)^2$ | $1.5171 \times 10^{-3}$ |
