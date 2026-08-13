# ============================================================
# Omega_parameter.py
# ============================================================
#
# Operational module for
#
#   A = Integral_0^1 F0(R) sinh(Psi(R)) R dR
#
#   B = Integral_0^1 F1(R) sinh(Psi(R)) R dR
#
#   Omega = A/(1-B)
#
# This module is intended for downstream use.
# It remains silent when imported.
#
# ============================================================

import numpy as np
from scipy.integrate import quad

import F0 as f0
import F1 as f1
import PoissonBoltzman as pb


# ============================================================
# NUMERICAL SETTINGS
# ============================================================

EPSABS_OMEGA = 1.0e-11
EPSREL_OMEGA = 1.0e-11
LIMIT_OMEGA = 300


# ============================================================
# INTEGRANDS
# ============================================================

def _integrand_A(R):
    return (
        f0.F0(R)
        * np.sinh(pb.Psi(R))
        * R
    )


def _integrand_B(R):
    return (
        f1.F1(R)
        * np.sinh(pb.Psi(R))
        * R
    )


# ============================================================
# COMPUTE A AND B
# ============================================================

def compute_A():
    """
    Compute

        A = Integral_0^1 F0(R) sinh(Psi(R)) R dR.
    """

    A, _ = quad(
        _integrand_A,
        0.0,
        1.0,
        epsabs=EPSABS_OMEGA,
        epsrel=EPSREL_OMEGA,
        limit=LIMIT_OMEGA
    )

    return A


def compute_B():
    """
    Compute

        B = Integral_0^1 F1(R) sinh(Psi(R)) R dR.
    """

    B, _ = quad(
        _integrand_B,
        0.0,
        1.0,
        epsabs=EPSABS_OMEGA,
        epsrel=EPSREL_OMEGA,
        limit=LIMIT_OMEGA
    )

    return B


# ============================================================
# MODULE-LEVEL VALUES
# ============================================================

A = compute_A()
B = compute_B()

denominator = 1.0 - B

if abs(denominator) < 1.0e-8:
    raise RuntimeError(
        "1 - B is numerically too small for a stable "
        "evaluation of Omega."
    )

Omega = A / denominator


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
    print("Omega parameter")
    print("---------------")
    print(f"A       = {A:.17e}")
    print(f"B       = {B:.17e}")
    print(f"1 - B   = {denominator:.17e}")
    print(f"Omega   = {Omega:.17e}")
