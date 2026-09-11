---
layout: default
title: Physical Model
math: true
---

## Physical Model and Governing Equations

The physical system consists of a cylindrical microchannel of radius $a$ and length $l$[cite: 3]. The inner surface of the microchannel is assumed to be negatively charged[cite: 3]. As the electrolyte flows through the channel, the mobile charge within the electrical double layer is transported in the axial direction, leading to charge accumulation at the channel ends and thereby generating a streaming potential[cite: 3].

<div style="text-align: center; margin: 30px 0;">
  <!-- Asegúrate de subir tu diagrama compilado como imagen a la carpeta de assets -->
  <img src="{{ '/assets/images/diagrama_microcanal.png' | relative_url }}" alt="Schematic representation of the cylindrical microchannel" style="max-width: 100%; border-radius: 8px; box-shadow: 0 4px 6px rgba(0,0,0,0.05);">
  <p style="color: #666; font-size: 0.9rem; margin-top: 10px;"><em>Figure 1: Schematic representation of the cylindrical microchannel.</em></p>
</div>

### Dimensionless Variables

To generalize the solution, the governing equations are transformed into a dimensionless form. The mathematical formulation uses the dimensionless variables defined in the following table:

| Dimensionless variable | Definition | Description |
| :--- | :--- | :--- |
| $Z$ | $z/l$ | Where $z$ is the axial coordinate and $l$ is the characteristic channel length.[cite: 3] |
| $R$ | $r/a$ | Where $r$ is the radial coordinate and $a$ is the radius of the microchannel.[cite: 3] |
| $V_Z$ | $v_z/J_w$ | Where $v_z$ is the axial velocity component and $J_w$ is defined via the volumetric flow rate as $Q=A_m J_w$.[cite: 3] |
| $V_R$ | $v_r/v_{r,c}$ | Where $v_r$ is the radial velocity component and $v_{r,c}$ is the characteristic radial velocity scale.[cite: 3] |
| $\Psi$ | $\psi/\zeta_T$ | Where $\psi$ is the electric double-layer (EDL) potential and $\zeta_T$ is the thermal potential.[cite: 3] |
| $\Phi$ | $\phi/\phi_c$ | Where $\phi$ is the streaming potential and $\phi_c$ is the characteristic streaming potential scale.[cite: 3] |
| $\widetilde{\Pi}$ | $\widetilde{p}/P_{\mathrm{OS,D,b}}$ | Where $\widetilde{p}$ is the modified pressure and $P_{\mathrm{OS,D,b}}$ is the bulk osmotic pressure reference scale.[cite: 3] |

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
