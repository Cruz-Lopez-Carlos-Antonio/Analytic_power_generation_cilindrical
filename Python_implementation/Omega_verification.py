# ============================================================
# Omega_verification.py
# ============================================================
#
# Validation of
#
#   A = Integral_0^1 F0(R) sinh(Psi(R)) R dR
#
#   B = Integral_0^1 F1(R) sinh(Psi(R)) R dR
#
#   Omega = A/(1-B)
#
# followed by the intrinsic consistency check
#
#   OmegaCheck =
#       Integral_0^1
#       [F0(R) + Omega F1(R)]
#       sinh(Psi(R)) R dR.
#
# The purpose of this script is validation, not operational use.
# A smaller silent module can be created afterwards.
#
# ============================================================


from scipy.integrate import quad

import F0 as f0
import F1 as f1
import PoissonBoltzman as pb


# ============================================================
# NUMERICAL SETTINGS
# ============================================================
#
# The whole Python chain works in double precision.
# The quadrature tolerances should therefore be demanding,
# but not interpreted as proof that the upstream BVP/ODE
# solutions contain the same number of reliable digits.
#
# ============================================================

EPSABS_OMEGA = 1.0e-11
EPSREL_OMEGA = 1.0e-11
LIMIT_OMEGA = 300


# ============================================================
# MATHEMATICA REFERENCE
# ============================================================
#
# Validated value from Omega_verification.wl.
# This is used ONLY for transversal comparison.
# It does not enter the computation of A, B, Omega or OmegaCheck.
#
# ============================================================

OMEGA_MATHEMATICA = 0.225503525392794229


# ============================================================
# INTEGRANDS
# ============================================================

def integrand_A(R):
    """
    Integrand for

        A = Integral_0^1 F0(R) sinh(Psi(R)) R dR.
    """

    import numpy as np

    return (
        f0.F0(R)
        * np.sinh(pb.Psi(R))
        * R
    )


def integrand_B(R):
    """
    Integrand for

        B = Integral_0^1 F1(R) sinh(Psi(R)) R dR.
    """

    import numpy as np

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
    Compute A and QUADPACK's absolute error estimate.
    """

    return quad(
        integrand_A,
        0.0,
        1.0,
        epsabs=EPSABS_OMEGA,
        epsrel=EPSREL_OMEGA,
        limit=LIMIT_OMEGA
    )


def compute_B():
    """
    Compute B and QUADPACK's absolute error estimate.
    """

    return quad(
        integrand_B,
        0.0,
        1.0,
        epsabs=EPSABS_OMEGA,
        epsrel=EPSREL_OMEGA,
        limit=LIMIT_OMEGA
    )


# ============================================================
# COMPUTE OMEGA
# ============================================================

def compute_Omega(A=None, B=None):
    """
    Compute

        Omega = A/(1-B).

    If A or B are not supplied, they are computed here.
    """

    if A is None:
        A, _ = compute_A()

    if B is None:
        B, _ = compute_B()

    denominator = 1.0 - B

    if abs(denominator) < 1.0e-8:
        raise RuntimeError(
            "1 - B is numerically too small for a stable "
            "evaluation of Omega."
        )

    Omega = A / denominator

    return Omega, denominator


# ============================================================
# RECONSTRUCT F FOR THE INTRINSIC CHECK
# ============================================================

def F_verification(R, Omega):
    """
    Reconstructed semianalytic function

        F(R) = F0(R) + Omega F1(R).
    """

    return (
        f0.F0(R)
        + Omega * f1.F1(R)
    )


def omega_check_integrand(R, Omega):
    """
    Integrand of the defining functional of Omega.
    """

    import numpy as np

    return (
        F_verification(R, Omega)
        * np.sinh(pb.Psi(R))
        * R
    )


def compute_OmegaCheck(Omega):
    """
    Compute

        OmegaCheck =
            Integral_0^1
            F(R) sinh(Psi(R)) R dR.
    """

    return quad(
        lambda R: omega_check_integrand(R, Omega),
        0.0,
        1.0,
        epsabs=EPSABS_OMEGA,
        epsrel=EPSREL_OMEGA,
        limit=LIMIT_OMEGA
    )


# ============================================================
# VALIDATION REPORT
# ============================================================

def validate_Omega():
    """
    Compute A, B, Omega and OmegaCheck, then compare:

      1. OmegaCheck vs Omega
         -> intrinsic consistency check;

      2. Omega(Python) vs Omega(Mathematica)
         -> transversal language-to-language comparison.
    """

    A, A_error = compute_A()
    B, B_error = compute_B()

    Omega, denominator = compute_Omega(
        A=A,
        B=B
    )

    OmegaCheck, OmegaCheck_error = compute_OmegaCheck(
        Omega
    )

    # Intrinsic comparison
    intrinsic_abs_error = abs(
        OmegaCheck - Omega
    )

    intrinsic_ape = (
        100.0
        * intrinsic_abs_error
        / abs(Omega)
        if Omega != 0.0
        else float("nan")
    )

    # Mathematica vs Python
    transversal_abs_error = abs(
        Omega - OMEGA_MATHEMATICA
    )

    transversal_ape = (
        100.0
        * transversal_abs_error
        / abs(OMEGA_MATHEMATICA)
        if OMEGA_MATHEMATICA != 0.0
        else float("nan")
    )

    print("")
    print("==============================================")
    print("A, B AND OMEGA VERIFICATION — PYTHON")
    print("==============================================")

    print("")
    print(f"A                  = {A:.17e}")
    print(f"quad error(A)      = {A_error:.6e}")

    print("")
    print(f"B                  = {B:.17e}")
    print(f"quad error(B)      = {B_error:.6e}")

    print("")
    print(f"1 - B              = {denominator:.17e}")

    print("")
    print(f"Omega              = {Omega:.17e}")
    print(f"OmegaCheck         = {OmegaCheck:.17e}")
    print(f"quad error(check)  = {OmegaCheck_error:.6e}")

    print("")
    print("Intrinsic consistency")
    print("---------------------")
    print(
        f"Absolute error     = "
        f"{intrinsic_abs_error:.6e}"
    )
    print(
        f"APE (%)            = "
        f"{intrinsic_ape:.6e}"
    )

    print("")
    print("Mathematica vs Python")
    print("---------------------")
    print(
        f"Omega Mathematica  = "
        f"{OMEGA_MATHEMATICA:.17e}"
    )
    print(
        f"Absolute error     = "
        f"{transversal_abs_error:.6e}"
    )
    print(
        f"APE (%)            = "
        f"{transversal_ape:.6e}"
    )

    print("")
    print("==============================================")

    return {
        "A": A,
        "A_quad_error": A_error,
        "B": B,
        "B_quad_error": B_error,
        "denominator": denominator,
        "Omega": Omega,
        "OmegaCheck": OmegaCheck,
        "OmegaCheck_quad_error": OmegaCheck_error,
        "intrinsic_abs_error": intrinsic_abs_error,
        "intrinsic_APE_percent": intrinsic_ape,
        "Omega_Mathematica": OMEGA_MATHEMATICA,
        "transversal_abs_error": transversal_abs_error,
        "transversal_APE_percent": transversal_ape,
    }


# ============================================================
# DIRECT EXECUTION
# ============================================================

if __name__ == "__main__":
    validate_Omega()
