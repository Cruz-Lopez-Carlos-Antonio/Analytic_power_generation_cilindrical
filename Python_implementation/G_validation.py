# ============================================================
# G_verification.py
# ============================================================
#
# Verification of the dimensionless modified-pressure gradient
#
#     G = 1 / [2 Integral_0^1 R F(R) dR].
#
# The script checks:
#
#   1. the flow integral IF;
#   2. non-degeneracy of the denominator;
#   3. the resulting pressure gradient G;
#   4. the flow-rate normalization;
#   5. agreement with Mathematica reference values.
#
# ============================================================

import numpy as np
from scipy.integrate import quad

import F as f


# ============================================================
# NUMERICAL SETTINGS
# ============================================================

EPSABS_G = 1.0e-12
EPSREL_G = 1.0e-12
QUAD_LIMIT_G = 300

DEGENERACY_TOL = 1.0e-12


# ============================================================
# MATHEMATICA REFERENCE VALUES
# ============================================================

IF_MATHEMATICA = -18.4404755931094882
G_MATHEMATICA = -0.0271142681475542412


# ============================================================
# FLOW INTEGRAND
# ============================================================

def flow_integrand(R):
    """
    Integrand for

        IF = Integral_0^1 R F(R) dR.
    """

    return R * f.F(R)


# ============================================================
# FLOW INTEGRAL
# ============================================================

IF, IF_quad_error = quad(
    flow_integrand,
    0.0,
    1.0,
    epsabs=EPSABS_G,
    epsrel=EPSREL_G,
    limit=QUAD_LIMIT_G
)


# ============================================================
# NON-DEGENERACY CHECK
# ============================================================

if not np.isfinite(IF):
    raise RuntimeError(
        "The flow integral IF is not finite."
    )

if abs(IF) < DEGENERACY_TOL:
    raise RuntimeError(
        "The denominator defining G is numerically degenerate."
    )


# ============================================================
# PRESSURE GRADIENT
# ============================================================

G = 1.0 / (2.0 * IF)


# ============================================================
# FLOW-RATE NORMALIZATION CHECK
# ============================================================

def normalized_flow_integrand(R):
    """
    Integrand for the normalized flow-rate condition

        Integral_0^1 R G F(R) dR = 1/2.
    """

    return R * G * f.F(R)


flow_check, flow_quad_error = quad(
    normalized_flow_integrand,
    0.0,
    1.0,
    epsabs=EPSABS_G,
    epsrel=EPSREL_G,
    limit=QUAD_LIMIT_G
)

flow_target = 0.5
flow_residual = flow_check - flow_target


# ============================================================
# PYTHON–MATHEMATICA COMPARISON
# ============================================================

IF_abs_error = abs(
    IF - IF_MATHEMATICA
)

G_abs_error = abs(
    G - G_MATHEMATICA
)

IF_APE = (
    IF_abs_error
    / abs(IF_MATHEMATICA)
    * 100.0
)

G_APE = (
    G_abs_error
    / abs(G_MATHEMATICA)
    * 100.0
)


# ============================================================
# OUTPUT
# ============================================================

print("")
print("============================================================")
print("Verification of dimensionless pressure gradient G")
print("============================================================")

print("")
print("Flow integral")
print("-------------")
print(f"IF              = {IF:.17e}")
print(f"quad error est. = {IF_quad_error:.3e}")

print("")
print("Pressure gradient")
print("-----------------")
print(f"G               = {G:.17e}")

print("")
print("Flow-rate normalization")
print("-----------------------")
print(f"Computed        = {flow_check:.17e}")
print(f"Target          = {flow_target:.17e}")
print(f"Abs. residual   = {abs(flow_residual):.3e}")
print(f"quad error est. = {flow_quad_error:.3e}")

print("")
print("Sign checks")
print("-----------")
print(f"Sign(IF)        = {np.sign(IF):.0f}")
print(f"Sign(G)         = {np.sign(G):.0f}")

print("")
print("Python vs Mathematica")
print("---------------------")

print(f"IF Mathematica  = {IF_MATHEMATICA:.17e}")
print(f"IF Python       = {IF:.17e}")
print(f"IF abs. error   = {IF_abs_error:.3e}")
print(f"IF APE (%)      = {IF_APE:.3e}")

print("")

print(f"G Mathematica   = {G_MATHEMATICA:.17e}")
print(f"G Python        = {G:.17e}")
print(f"G abs. error    = {G_abs_error:.3e}")
print(f"G APE (%)       = {G_APE:.3e}")

print("")
print("============================================================")
