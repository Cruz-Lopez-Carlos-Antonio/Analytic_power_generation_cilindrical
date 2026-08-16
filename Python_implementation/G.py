# ============================================================
# G.py
# ============================================================
#
# Operational module for the dimensionless modified-pressure
# gradient
#
#     G = 1 / [2 Integral_0^1 R F(R) dR].
#
# This module is intended for downstream use.
# It remains silent when imported.
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


# ============================================================
# FLOW INTEGRAND
# ============================================================

def _flow_integrand(R):
    """
    Integrand for

        IF = Integral_0^1 R F(R) dR.
    """

    return R * f.F(R)


# ============================================================
# FLOW INTEGRAL
# ============================================================

IF, _IF_quad_error = quad(
    _flow_integrand,
    0.0,
    1.0,
    epsabs=EPSABS_G,
    epsrel=EPSREL_G,
    limit=QUAD_LIMIT_G
)


# ============================================================
# DIMENSIONLESS PRESSURE GRADIENT
# ============================================================

G = 1.0 / (2.0 * IF)


# ============================================================
# DIRECT EXECUTION
# ============================================================
#
# Importing this module produces no output.
# Running it directly performs only a minimal operational check.
#
# ============================================================

if __name__ == "__main__":

    print("")
    print("G operational check")
    print("-------------------")

    print(f"IF = {IF:.17e}")
    print(f"G  = {G:.17e}")
