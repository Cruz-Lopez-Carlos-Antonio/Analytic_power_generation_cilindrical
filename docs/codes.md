---
layout: default
title: Codes
math: true
---

# Semianalytical Parametric Solvers

The computational core of this repository focuses on solving two primary physical variables of interest: the average streaming potential, $\Phi_{\mathrm{av}}$, and the conversion efficiency, $\eta$, whose expressions are explicitly given in [the Semi-analytic Solution section](https://github.com/Cruz-Lopez-Carlos-Antonio/Analytic_power_generation_cilindrical/blob/main/docs/semianalytic.md#streaming-potential-and-conversion-efficiency) of this repository. These variables are evaluated as parametric functions, meaning their outputs depend directly on a defined set of physical and geometric inputs.. These variables are evaluated as parametric functions, meaning their outputs depend directly on a defined set of physical and geometric inputs. 

## 1. Average Streaming potential, $\Phi_{\mathrm{av}}$:

Depending on the specific numerical experiment, the general code automatically generates arrays of evaluations by treating different parameters as independent variables. Specifically, the solver computes the following three functional dependencies:

<div style="background:#f7f7f7; padding:15px; border-left:4px solid #4a90e2; border-radius:6px; margin:20px 0;">
$$
\begin{aligned}
\Phi_{\mathrm{av}} &= \Phi_{\mathrm{av}}(\delta, \Psi_s, \Gamma_1), \\
\Phi_{\mathrm{av}} &= \Phi_{\mathrm{av}}(\delta, \omega, \Gamma_2), \\
\Phi_{\mathrm{av}} &= \Phi_{\mathrm{av}}(\delta, \Lambda, \Gamma_3),
\end{aligned}
$$
</div>

where $\delta$ is the electrokinetic radius, $\Psi_s$ is the surface potential, $\omega$ is the viscoelectric parameter, and $\Lambda$ is the electrokinetic coupling parameter. The term $\Gamma_i$ represents the complementary set of fixed physical variables supplied by the base configuration (`Parameters.wl`), explicitly excluding the parameter currently being varied.

*   $\Gamma_1 = \{\omega, \Lambda, \Pi_D\}$
*   $\Gamma_2 = \{\Psi_s, \Lambda, \Pi_D\}$
*   $\Gamma_3 = \{\Psi_s, \omega, \Pi_D\}$

*(Note: Although the longitudinal pressure measure $\Pi_D$ is formally part of the parametric solver inputs, the average streaming potential is analytically independent of this parameter. Therefore, parametric sweeps involving $\Pi_D$ are programmatically skipped by the main script to save computational resources).*

### 1.1 Evaluated Intervals

By default, the automated sweeps are defined over the following numerical domains:

*   **Electrokinetic radius ($\delta$):** From $1$ to $30$ with step increments of $1$.
*   **Surface potential ($\Psi_s$):** $\{-1.3, -1.23249, -1.2, -1.1, -1.0, -0.9, -0.8, -0.7, -0.6, -0.5\}$.
*   **Viscoelectric parameter ($\omega$):** From $0$ to $0.005$ with increments of $0.0005$.
*   **Coupling parameter ($\Lambda$):** From $0.15$ to $1.05$ with increments of $0.15$.

<div style="border-left: 4px solid #f39c12; padding: 0.7em 1em; background: #fff7e6; margin: 20px 0;">
<b style="color:#c0392b;">⚠️ Important regarding numerical precision:</b><br>
When evaluating individual points using the solver module, any explicit parameter value supplied must be an <strong>exact number</strong> (integer or rational fraction, e.g., <code>45/100</code>). Do not supply machine-precision decimals (e.g., <code>0.45</code>), as the solver deliberately preserves infinite precision strings and will fail to execute inexact inputs.
</div>

### 1.2 Main Script and Dependencies

The automated parametric studies for $\Phi_{\mathrm{av}}$ are generated using the main script, available here:

👉 <a href="https://github.com/Cruz-Lopez-Carlos-Antonio/Analytic_power_generation_cilindrical/blob/main/Mathematica/Parametric/GeneratePhiAvParametricStudy.wl" target="_blank" rel="noopener noreferrer">GeneratePhiAvParametricStudy.wl</a>

To execute successfully, this main script relies on a hierarchical structure of dependencies. It primarily calls the central wrapper (`SemianalyticalParametricSolver_exact_inputs.wl`), which in turn loads the base parameters and all the validated analytical modules. 

<div style="background:#f1f7ff; padding:15px; border-left:4px solid #4a90e2; border-radius:8px; margin-top:20px;">
  <strong>📂 Required Dependency Modules:</strong>
  <ul style="margin-top: 10px; margin-bottom: 0;">
    <li><code>SemianalyticalParametricSolver_exact_inputs.wl</code> (Central Wrapper)</li>
    <li><code>Parameters.wl</code> (Base variables)</li>
    <li><code>PoissonBoltzmann_parametric.wl</code></li>
    <li><code>F_cc_parametric.wl</code></li>
    <li><code>Lambda_parameter_parametric.wl</code></li>
    <li><code>MR_parametric.wl</code></li>
    <li><code>HR_parametric.wl</code></li>
    <li><code>F0_parametric.wl</code></li>
    <li><code>F1_parametric.wl</code></li>
    <li><code>Omega_parameter_parametric.wl</code></li>
    <li><code>F_parametric.wl</code></li>
    <li><code>G_parametric.wl</code></li>
  </ul>
</div>

### 1.3 Modification of the Input Intervals

The numerical intervals used for the parametric studies are defined manually at the very end of the `GeneratePhiAvParametricStudy.wl` script. If user needs to evaluate a different range of physical parameters, you can directly modify these exact arrays within the code.

**Interval for $\Psi_s$:**
```mathematica
PhiAvPsiSValues = {-13/10,-123249/100000,
   -6/5,-11/10,-1,-9/10,-4/5,
   -7/10,-3/5,-1/2};
```

**Interval for $\Lambda$:**
```mathematica
PhiAvLambdaValues =
   Range[15, 105, 15]/100;
```
**Interval for $\omega$:**
```mathematica
PhiAvLambdaValues =
   Range[15, 105, 15]/100;
```





