# ============================================================
# F.py
# ============================================================
#
# Operational module for the reconstructed semianalytic function
#
#     F(R) = F0(R) + Omega F1(R)
#
# and its radial derivative
#
#     F'(R)
#       =
#     F0'(R)
#       +
#     Omega Lambda delta^2 Q(R).
#
# This module is intended for downstream use.
# It remains silent when imported.
#
# ============================================================

import numpy as np

import Parameters as par
import Lambda_parameter as lam
import F0 as f0
import F1 as f1
import Omega_parameter as om


# ============================================================
# PUBLIC FUNCTION: F(R)
# ============================================================

def F(R):
    """
    Reconstructed semianalytic function

        F(R) = F0(R) + Omega F1(R).

    Parameters
    ----------
    R : float or array_like
        Dimensionless radial coordinate.

        Required domain:
            0 <= R <= 1

    Returns
    -------
    float or ndarray
        F(R).
    """

    R_array = np.asarray(
        R,
        dtype=float
    )

    if np.any(
        (R_array < 0.0)
        | (R_array > 1.0)
    ):
        raise ValueError(
            "R must satisfy 0 <= R <= 1."
        )

    values = (
        f0.F0(R_array)
        + om.Omega * f1.F1(R_array)
    )

    # Enforce exact wall value:
    #
    #   F(1) = 0.
    #
    values = np.where(
        R_array == 1.0,
        0.0,
        values
    )

    if np.ndim(R) == 0:
        return float(values)

    return values


# ============================================================
# PUBLIC FUNCTION: F'(R)
# ============================================================

def FPrime(R):
    """
    Radial derivative of F(R):

        F'(R)
          =
        F0'(R)
          +
        Omega F1'(R),

    where

        F1'(R)
          =
        Lambda delta^2 Q(R).

    Regularity implies

        F'(0) = 0.
    """

    R_array = np.asarray(
        R,
        dtype=float
    )

    if np.any(
        (R_array < 0.0)
        | (R_array > 1.0)
    ):
        raise ValueError(
            "R must satisfy 0 <= R <= 1."
        )

    values = (
        f0.F0Prime(R_array)
        +
        om.Omega
        * lam.Lambda
        * par.delta**2
        * f1.QF1(R_array)
    )

    # Enforce exact symmetry-axis value:
    #
    #   F'(0) = 0.
    #
    values = np.where(
        R_array == 0.0,
        0.0,
        values
    )

    if np.ndim(R) == 0:
        return float(values)

    return values


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
    print("F(R) operational check")
    print("----------------------")

    print(f"Omega   = {om.Omega:.17e}")

    print("")
    print(f"F(0)    = {F(0.0):.15e}")
    print(f"F(0.5)  = {F(0.5):.15e}")
    print(f"F(1)    = {F(1.0):.15e}")

    print("")
    print(f"F'(0)   = {FPrime(0.0):.15e}")
