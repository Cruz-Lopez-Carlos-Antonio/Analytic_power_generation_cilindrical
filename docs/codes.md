---
layout: default
title: Codes
math: true
---

# Semianalytical Parametric Solvers

The computational core of this repository focuses on solving two primary physical variables of interest: the average streaming potential, $\Phi_{\mathrm{av}}$, and the conversion efficiency, $\eta$, whose expressions are explicitly given in Section [] of this repository. These variables are evaluated as parametric functions, meaning their outputs depend directly on a defined set of physical and geometric inputs. 

## 1. Average Streaming potential, $\Phi_{\mathrm{av}}$:

<div style="background:#f7f7f7; padding:15px; border-left:4px solid #4a90e2; border-radius:6px; margin:20px 0;">
$$
\Phi_{\mathrm{av}} = \Phi_{\mathrm{av}}(\delta, \Psi_s, \text{set})
$$
</div>

where $\delta$ is the electrokinetic radius, $\Psi_s$ is the surface potential, and $\text{set}$ represents a grouped array of additional input parameters (such as the viscoelectric parameter $\omega$, the coupling parameter $\Lambda$, and the longitudinal pressure measure $\Pi_D$) evaluated seamlessly by the solver. The corresponding efficiency solver $\eta(\delta, \Psi_s, \text{set})$ will be incorporated once its definitive module is validated.

### Main Script and Dependencies

The automated parametric studies for $\Phi_{\mathrm{av}}$ are generated using the main script, available here:

👉 <a href="https://github.com/Cruz-Lopez-Carlos-Antonio/Analytic_power_generation_cilindrical/blob/main/Mathematica/Parametric/GeneratePhiAvParametricStudy.wl" target="_blank" rel="noopener noreferrer">GeneratePhiAvParametricStudy.wl</a>

To execute successfully, this main script relies on a hierarchical structure of dependencies. It primarily calls the central wrapper (`SemianalyticalParametricSolver_exact_inputs.wl`), which in turn loads the base parameters and all the validated analytical modules. 

<div style="background:#f1f7ff; padding:15px; border-left:4px solid #4a90e2; border-radius:8px; margin-top:20px;">
  <strong>📂 Required Dependency Modules:</strong>
  <ul style="margin-top: 10px; margin-bottom: 0;">
    <li><code>SemianalyticalParametricSolver_exact_inputs.wl</code> (Central Wrapper)</li>
    <li><code>Parameters.wl</code> (Base variables)</li>
    <li><code>PoissonBoltzmann_parametric.wl</code></li>
    <li><code>F_cc_parametric.wl</code>[cite: 9]</li>
    <li><code>Lambda_parameter_parametric.wl</code></li>
    <li><code>MR_parametric.wl</code></li>
    <li><code>HR_parametric.wl</code></li>
    <li><code>F0_parametric.wl</code></li>
    <li><code>F1_parametric.wl</code></li>
    <li><code>Omega_parameter_parametric.wl</code>[cite: 9]</li>
    <li><code>F_parametric.wl</code>[cite: 9]</li>
    <li><code>G_parametric.wl</code>[cite: 9]</li>
  </ul>
</div>
