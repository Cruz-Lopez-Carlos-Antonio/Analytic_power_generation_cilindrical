---
layout: default
title: Overview
math: true
---
# Overview

The present repository documents the computational implementation and parametric solvers for the research:

<div style="display:inline-block; padding:0.6em 0.9em; margin:0.4em 0; background-color:#f5f5f5; border:1px solid #dddddd; border-radius:6px;" markdown="1">
***Power generation in cylindrical microchannels with high surface zeta potential: a survey on viscoelectric effects over streaming potentials and efficiency***
</div>

Developed by Sánchez Lozano, G., Cruz-López C.-A., and F. Méndez (2026).

## Wolfram Mathematica Codes

<div style="margin: 1.5rem 0; padding: 1.25rem 1.5rem; background:#fafafa; border-radius:10px; border:1px solid #e0e0e0;">
  <h3 style="margin-top:0; margin-bottom:0.75rem; font-size:1.1rem;">
    Implementations included in this repository:
  </h3>
  <ul style="margin:0; padding-left:1.2rem; list-style-type:disc;">
    <li style="margin:0.4rem 0;">
      <strong>Semianalytical Parametric Solver</strong><br/>
      <span>A general wrapper for evaluating physical parameters seamlessly.</span>
    </li>
    <li style="margin:0.4rem 0;">
      <strong>Direct Augmented Solver</strong><br/>
      <span>High-precision discrete evaluation of the system constraints.</span>
    </li>
    <li style="margin:0.4rem 0;">
      <strong>Efficiency ($\eta$) and Potential ($\Phi_{av}$) Studies</strong><br/>
      <span>Automated parametric sweeps for numerical validation.</span>
    </li>
  </ul>
</div>

<div style="text-align: center; margin: 30px 0;">
  <!-- Sube la imagen compilada de tu esquema de dependencias de Wolfram a la carpeta assets/images/ -->
  <img src="{{ '/assets/images/wolfram_solver_workflow.png' | relative_url }}" alt="Hierarchical dependency scheme of the parameterized semianalytical solver" style="max-width: 100%; border-radius: 8px; box-shadow: 0 4px 6px rgba(0,0,0,0.05);">
  <p style="color: #666; font-size: 0.9rem; margin-top: 10px;"><em>Figure: Hierarchical workflow of the semianalytical parametric modules.</em></p>
</div>

# Map of the site 
You can use the navigation bar above, or the following buttons, to explore:

<div style="display:grid; grid-template-columns:repeat(auto-fit,minmax(220px,1fr)); gap:0.9rem; margin:1rem 0;">

  <a href="{{ '/equations.html' | relative_url }}" style="text-decoration:none; color:inherit;">
    <div style="background:#fafafa; border:1px solid #e0e0e0; border-radius:10px; padding:0.75rem 0.9rem; cursor:pointer;">
      <strong>Physical description and<br>governing equations</strong><br/>
      <span style="font-size:0.95rem; color:#555;">
        Main analytical expressions and key formulae used in the article.
      </span>
    </div>
  </a>

  <a href="{{ '/codes.html' | relative_url }}" style="text-decoration:none; color:inherit;">
    <div style="background:#fafafa; border:1px solid #e0e0e0; border-radius:10px; padding:0.75rem 0.9rem; cursor:pointer;">
      <strong>Codes</strong><br/>
      <span style="font-size:0.95rem; color:#555;">
        Summary of the Wolfram Mathematica scripts, interfaces, and numerical settings.
      </span>
    </div>
  </a>

  <a href="{{ '/validation.html' | relative_url }}" style="text-decoration:none; color:inherit;">
    <div style="background:#fafafa; border:1px solid #e0e0e0; border-radius:10px; padding:0.75rem 0.9rem; cursor:pointer;">
      <strong>Validation</strong><br/>
      <span style="font-size:0.95rem; color:#555;">
        Comparison against the numerical benchmarks and parametric studies.
      </span>
    </div>
  </a>

  <a href="{{ '/about.html' | relative_url }}" style="text-decoration:none; color:inherit;">
    <div style="background:#fafafa; border:1px solid #e0e0e0; border-radius:10px; padding:0.75rem 0.9rem; cursor:pointer;">
      <strong>About</strong><br/>
      <span style="font-size:0.95rem; color:#555;">
        Authorship, affiliations, and financial support acknowledgements.
      </span>
    </div>
  </a>

</div>

<div style="border-left: 4px solid #f39c12; padding: 0.7em 1em; background: #fff7e6; margin-top: 25px;">
<b style="color:#c0392b;">⚠️  Important:</b><br>
For more details on the derivation of the physical model, its computational implementation, as well as a comprehensive numerical analysis, please see the manuscript mentioned above.
</div>
