# ============================================================
# F_verification.py
# ============================================================
#
# Validation of the reconstructed semianalytic function
#
#     F(R) = F0(R) + Omega F1(R)
#
# using the already validated operational modules:
#
#     F0.py
#     F1.py
#     Omega_parameter.py
#
# The script checks:
#
#   1. F(1) = 0
#   2. F'(0) = 0
#   3. Omega = Integral_0^1 F(R) sinh(Psi(R)) R dR
#   4. residual of the integrated Eq. (50):
#
#        M(R) R F'(R)
#          =
#        Pi_D R^2/2
#          +
#        Omega Lambda delta^2 H(R)
#
# ============================================================

import numpy as np
from scipy.integrate import quad

import Parameters as par
import PoissonBoltzman as pb
import MR as mr
import HR as hr
import Lambda_parameter as lam
import F0 as f0
import F1 as f1
import Omega_parameter as om


# ============================================================
# NUMERICAL SETTINGS
# ============================================================

EPSABS_F = 1.0e-11
EPSREL_F = 1.0e-11
LIMIT_F = 300

N_RESIDUAL = 1001


# ============================================================
# RECONSTRUCTED F(R)
# ============================================================

def F(R):
    """
    Reconstructed semianalytic function

        F(R) = F0(R) + Omega F1(R).
    """

    return (
        f0.F0(R)
        + om.Omega * f1.F1(R)
    )


# ============================================================
# DERIVATIVE F'(R)
# ============================================================

def FPrime(R):
    """
    Derivative of F(R):

        F'(R)
          =
        F0'(R)
          +
        Omega F1'(R),

    with

        F0'(R)
          =
        Pi_D R/(2 M(R)),

    and

        F1'(R)
          =
        Lambda delta^2 Q(R).
    """

    R_array = np.asarray(
        R,
        dtype=float
    )

    values = (
        f0.F0Prime(R_array)
        +
        om.Omega
        * lam.Lambda
        * par.delta**2
        * f1.QF1(R_array)
    )

    if np.ndim(R) == 0:
        return float(values)

    return values


# ============================================================
# OMEGA FUNCTIONAL CHECK
# ============================================================

def omega_functional_integrand(R):
    """
    Integrand of the defining functional

        Omega =
        Integral_0^1 F(R) sinh(Psi(R)) R dR.
    """

    return (
        F(R)
        * np.sinh(pb.Psi(R))
        * R
    )


def compute_omega_check():
    """
    Recompute Omega directly from its defining integral.
    """

    return quad(
        omega_functional_integrand,
        0.0,
        1.0,
        epsabs=EPSABS_F,
        epsrel=EPSREL_F,
        limit=LIMIT_F
    )


# ============================================================
# RESIDUAL OF THE INTEGRATED EQUATION
# ============================================================

def integrated_residual(R):
    """
    Residual of

        M(R) R F'(R)
          =
        Pi_D R^2/2
          +
        Omega Lambda delta^2 H(R).

    This form avoids numerical differentiation of the flux.
    """

    R_array = np.asarray(
        R,
        dtype=float
    )

    lhs = (
        mr.MR(R_array)
        * R_array
        * FPrime(R_array)
    )

    rhs = (
        f0.PiD
        * R_array**2
        / 2.0
        +
        om.Omega
        * lam.Lambda
        * par.delta**2
        * hr.HR(R_array)
    )

    residual = lhs - rhs

    if np.ndim(R) == 0:
        return float(residual)

    return residual


# ============================================================
# VALIDATION REPORT
# ============================================================

def validate_F():
    """
    Perform the main checks for F(R).
    """

    # Boundary and symmetry conditions
    F_wall = F(1.0)
    Fprime_axis = FPrime(0.0)

    # Omega functional
    OmegaCheck, OmegaCheck_quad_error = compute_omega_check()

    Omega_abs_error = abs(
        OmegaCheck - om.Omega
    )

    Omega_ape = (
        100.0
        * Omega_abs_error
        / abs(om.Omega)
        if om.Omega != 0.0
        else float("nan")
    )

    # Integrated-equation residual on a dense grid
    R_grid = np.linspace(
        0.0,
        1.0,
        N_RESIDUAL
    )

    residual = integrated_residual(
        R_grid
    )

    residual_linf = np.max(
        np.abs(residual)
    )

    residual_l2 = np.sqrt(
        np.mean(residual**2)
    )

    # Selected values
    R_sample = np.array([
        0.0,
        0.25,
        0.50,
        0.75,
        1.0
    ])

    F_sample = F(
        R_sample
    )

    print("")
    print("==============================================")
    print("F(R) VERIFICATION — PYTHON")
    print("==============================================")

    print("")
    print(f"Omega              = {om.Omega:.17e}")

    print("")
    print("Boundary / symmetry conditions")
    print("------------------------------")
    print(f"F(1)               = {F_wall:.17e}")
    print(f"F'(0)              = {Fprime_axis:.17e}")

    print("")
    print("Omega functional check")
    print("----------------------")
    print(f"OmegaCheck         = {OmegaCheck:.17e}")
    print(f"quad error(check)  = {OmegaCheck_quad_error:.6e}")
    print(f"Absolute error     = {Omega_abs_error:.6e}")
    print(f"APE (%)            = {Omega_ape:.6e}")

    print("")
    print("Integrated Eq. (50) residual")
    print("----------------------------")
    print(f"L_inf residual     = {residual_linf:.6e}")
    print(f"L2 residual        = {residual_l2:.6e}")

    print("")
    print("Selected F(R) values")
    print("--------------------")
    print("       R                   F(R)")
    print("----------------------------------------")

    for R_i, F_i in zip(
        R_sample,
        F_sample
    ):
        print(
            f"{R_i:8.2f}     "
            f"{F_i: .15e}"
        )

    print("")
    print("==============================================")

    return {
        "Omega": om.Omega,
        "F_wall": F_wall,
        "Fprime_axis": Fprime_axis,
        "OmegaCheck": OmegaCheck,
        "OmegaCheck_quad_error": OmegaCheck_quad_error,
        "Omega_abs_error": Omega_abs_error,
        "Omega_APE_percent": Omega_ape,
        "residual_Linf": residual_linf,
        "residual_L2": residual_l2,
        "R_sample": R_sample,
        "F_sample": F_sample,
    }


# ============================================================
# DIRECT EXECUTION
# ============================================================

if __name__ == "__main__":
    validate_F()
