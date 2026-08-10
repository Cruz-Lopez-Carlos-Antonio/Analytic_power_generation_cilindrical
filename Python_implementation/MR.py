import numpy as np

import Parameters as par
import PoissonBoltzman as pb


# ============================================================
# M(R) — VISCOELECTRIC FACTOR
# ============================================================
#
# We define
#
#   M(R) = exp[omega * (Psi'(R))^2]
#
# where:
#
#   omega
#
# is obtained from Parameters.py, and
#
#   Psi'(R)
#
# is obtained from the validated Poisson-Boltzmann solution.
#
# By symmetry,
#
#   Psi'(0) = 0,
#
# therefore
#
#   M(0) = 1.
#
# ============================================================


def MR(R):
    """
    Dimensionless viscoelectric factor

        M(R) = exp[omega * (Psi'(R))^2].

    Parameters
    ----------
    R : float or array_like
        Dimensionless radial coordinate.

        Required domain:
            0 <= R <= 1

    Returns
    -------
    float or ndarray
        M(R).
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

    psi_prime = np.asarray(
        pb.PsiPrime(R_array),
        dtype=float
    )

    M_values = np.exp(
        par.omega * psi_prime**2
    )

    # Enforce the exact symmetry-axis value
    # for numerical consistency:
    #
    #   Psi'(0) = 0  ->  M(0) = 1.
    #
    M_values = np.where(
        R_array == 0.0,
        1.0,
        M_values
    )

    if np.ndim(R) == 0:
        return float(M_values)

    return M_values


# ============================================================
# VALIDATION VALUES
# ============================================================

def sample_MR():
    """
    Evaluate M(R) at 10 uniformly spaced radial positions
    for cross-validation with Mathematica and collaborator data.
    """

    R_sample = np.linspace(
        0.0,
        1.0,
        10
    )

    M_sample = MR(
        R_sample
    )

    print(
        "\nSelected M(R) values"
    )

    print(
        "------------------------------------------------"
    )

    print(
        "          R                         M(R)"
    )

    print(
        "------------------------------------------------"
    )

    for R_i, M_i in zip(
        R_sample,
        M_sample
    ):
        print(
            f"{R_i: .15f}     "
            f"{M_i: .17e}"
        )

    return (
        R_sample,
        M_sample
    )


# ============================================================
# BASIC VALIDATION
# ============================================================

def validate_MR():
    """
    Basic numerical checks for M(R).
    """

    R_test = np.linspace(
        0.0,
        1.0,
        1001
    )

    M_test = MR(
        R_test
    )

    print(
        "\nM(R) validation"
    )

    print(
        "----------------"
    )

    print(
        f"M(0)          = {MR(0.0):.17e}"
    )

    print(
        f"M(1)          = {MR(1.0):.17e}"
    )

    print(
        f"min M(R)      = {np.min(M_test):.17e}"
    )

    print(
        f"M(R) > 0      = {np.all(M_test > 0.0)}"
    )


# ============================================================
# DIRECT EXECUTION
# ============================================================
#
# The validation output is produced only when this file
# is executed directly.
#
# Importing
#
#   import MR
#
# remains silent.
#
# ============================================================

if __name__ == "__main__":

    validate_MR()

    sample_MR()
