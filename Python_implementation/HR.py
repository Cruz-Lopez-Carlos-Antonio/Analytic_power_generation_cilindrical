import numpy as np

from scipy.integrate import solve_ivp

import PoissonBoltzman as pb


# ============================================================
# H(R) — AUXILIARY ODE
# ============================================================
#
# We define
#
#   H(R) = integral_0^R tau*sinh(Psi(tau)) dtau
#
# therefore
#
#   H'(R) = R*sinh(Psi(R)),
#
# with
#
#   H(0) = 0.
#
# The auxiliary ODE is solved once and the resulting
# dense solution is cached for subsequent evaluations.
#
# This avoids nested quadratures in later calculations.
#
# ============================================================


# ============================================================
# AUXILIARY ODE
# ============================================================

def h_ode_system(R, y):
    """
    First-order auxiliary ODE for H(R):

        H'(R) = R*sinh(Psi(R)).
    """

    return [
        R * np.sinh(
            pb.Psi(R)
        )
    ]


# ============================================================
# CACHED ODE SOLUTION
# ============================================================

_h_solution = None


def get_h_solution():
    """
    Return the cached auxiliary-ODE solution.

    The ODE is solved only on the first call.
    """

    global _h_solution

    if _h_solution is None:

        _h_solution = solve_ivp(
            h_ode_system,
            (0.0, 1.0),
            [0.0],
            method="DOP853",
            rtol=1.0e-11,
            atol=1.0e-13,
            dense_output=True
        )

        if not _h_solution.success:
            raise RuntimeError(
                "H(R) auxiliary ODE solver failed: "
                + _h_solution.message
            )

    return _h_solution


# ============================================================
# PUBLIC FUNCTION: H(R)
# ============================================================

def HR(R):
    """
    Dimensionless auxiliary function

        H(R) = integral_0^R tau*sinh(Psi(tau)) dtau.

    Parameters
    ----------
    R : float or array_like
        Dimensionless radial coordinate.

        Required domain:
            0 <= R <= 1

    Returns
    -------
    float or ndarray
        H(R).
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

    solution = get_h_solution()

    H_values = solution.sol(
        R_array
    )[0]

    # Enforce the exact axis value
    #
    #   H(0) = 0.
    #
    H_values = np.where(
        R_array == 0.0,
        0.0,
        H_values
    )

    if np.ndim(R) == 0:
        return float(H_values)

    return H_values


# ============================================================
# DIRECT EXECUTION
# ============================================================
#
# Nothing is printed when this module is imported.
#
# Running this file directly performs only a minimal
# operational check.
#
# ============================================================

if __name__ == "__main__":

    print(
        "H(R) operational check"
    )

    print(
        "----------------------"
    )

    print(
        f"H(0)   = {HR(0.0):.15e}"
    )

    print(
        f"H(0.5) = {HR(0.5):.15e}"
    )

    print(
        f"H(1)   = {HR(1.0):.15e}"
    )
