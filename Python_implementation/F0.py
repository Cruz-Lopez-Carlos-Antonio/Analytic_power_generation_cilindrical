import numpy as np

from scipy.integrate import solve_ivp

import Parameters as par
import MR as mr


# ============================================================
# F0(R) — OPERATIONAL IMPLEMENTATION
# ============================================================
#
# We define
#
#   F0(R) =
#       -Integral_R^1 [PiD s/(2 M(s))] ds
#
# where
#
#   PiD = alpha xi^2.
#
# Instead of evaluating the quadrature repeatedly,
# we solve the equivalent auxiliary ODE
#
#   F0'(R) = PiD R/(2 M(R)),
#
# with
#
#   F0(1) = 0.
#
# The ODE is solved once and the dense solution is cached.
#
# ============================================================


# ============================================================
# PARAMETER
# ============================================================

PiD = (
    par.alpha
    * par.xi**2
)


# ============================================================
# NUMERICAL SETTINGS
# ============================================================

RTOL_F0 = 1.0e-11
ATOL_F0 = 1.0e-13


# ============================================================
# AUXILIARY ODE
# ============================================================

def f0_ode_system(R, y):
    """
    Auxiliary ODE for F0(R):

        F0'(R) = PiD R/(2 M(R)).
    """

    return [
        PiD
        * R
        / (
            2.0
            * mr.MR(R)
        )
    ]


# ============================================================
# CACHED ODE SOLUTION
# ============================================================

_f0_solution = None


def get_f0_solution():
    """
    Return the cached auxiliary-ODE solution for F0(R).

    Since

        F0(1) = 0,

    solve_ivp integrates backward from R = 1 to R = 0.

    The ODE is solved only on the first call.
    """

    global _f0_solution

    if _f0_solution is None:

        _f0_solution = solve_ivp(
            f0_ode_system,
            (1.0, 0.0),
            [0.0],
            method="DOP853",
            rtol=RTOL_F0,
            atol=ATOL_F0,
            dense_output=True
        )

        if not _f0_solution.success:
            raise RuntimeError(
                "F0(R) auxiliary ODE solver failed: "
                + _f0_solution.message
            )

    return _f0_solution


# ============================================================
# PUBLIC FUNCTION: F0(R)
# ============================================================

def F0(R):
    """
    Dimensionless auxiliary function F0(R).

    Parameters
    ----------
    R : float or array_like
        Dimensionless radial coordinate.

        Required domain:
            0 <= R <= 1

    Returns
    -------
    float or ndarray
        F0(R).
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

    solution = get_f0_solution()

    F0_values = solution.sol(
        R_array
    )[0]

    # Enforce exact boundary condition:
    #
    #   F0(1) = 0.
    #
    F0_values = np.where(
        R_array == 1.0,
        0.0,
        F0_values
    )

    if np.ndim(R) == 0:
        return float(F0_values)

    return F0_values


# ============================================================
# PUBLIC FUNCTION: F0'(R)
# ============================================================

def F0Prime(R):
    """
    Radial derivative of F0(R):

        F0'(R) = PiD R/(2 M(R)).
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
        PiD
        * R_array
        / (
            2.0
            * mr.MR(R_array)
        )
    )

    if np.ndim(R) == 0:
        return float(values)

    return values


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

    print("")
    print(
        "F0(R) operational check"
    )

    print(
        "-----------------------"
    )

    print(
        f"PiD       = {PiD:.15e}"
    )

    print(
        f"F0(0)     = {F0(0.0):.15e}"
    )

    print(
        f"F0(0.5)   = {F0(0.5):.15e}"
    )

    print(
        f"F0(1)     = {F0(1.0):.15e}"
    )

    print(
        f"F0'(0)    = {F0Prime(0.0):.15e}"
    )
