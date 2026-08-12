import numpy as np

from scipy.integrate import solve_ivp

import Parameters as par
import Lambda_parameter as lam
import PoissonBoltzman as pb
import HR as hr
import MR as mr


# ============================================================
# F1(R) — OPERATIONAL IMPLEMENTATION
# ============================================================
#
# We define
#
#   F1(R) =
#       -Lambda delta^2 Integral_R^1 Q(s) ds
#
# with
#
#   Q(R) = H(R)/(R M(R)).
#
# Instead of evaluating the quadrature repeatedly,
# we solve the equivalent auxiliary ODE:
#
#   F1'(R) = Lambda delta^2 Q(R),
#
# with boundary condition
#
#   F1(1) = 0.
#
# The ODE is solved once and the dense solution is cached.
#
# ============================================================


# ============================================================
# NUMERICAL SETTINGS
# ============================================================

EPS_Q = 1.0e-8

RTOL_F1 = 1.0e-11
ATOL_F1 = 1.0e-13


# ============================================================
# AUXILIARY FUNCTION Q(R)
# ============================================================

def QF1(R):
    """
    Auxiliary function

        Q(R) = H(R)/(R M(R)).

    The apparent singularity at R = 0 is removable.

    Near the symmetry axis,

        Q(R)
        = [sinh(Psi(0))/2] R + O(R^3).

    Parameters
    ----------
    R : float or array_like
        Dimensionless radial coordinate.

        Required domain:
            0 <= R <= 1

    Returns
    -------
    float or ndarray
        Q(R).
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

    Q_values = np.empty_like(
        R_array,
        dtype=float
    )

    mask_axis = (
        R_array == 0.0
    )

    mask_near = (
        (R_array > 0.0)
        & (R_array < EPS_Q)
    )

    mask_regular = (
        R_array >= EPS_Q
    )

    # Exact removable value:
    #
    #   Q(0) = 0.
    #
    Q_values[mask_axis] = 0.0

    # Local regular expansion:
    #
    #   Q(R) ~ 0.5*sinh(Psi(0))*R.
    #
    if np.any(mask_near):

        psi_center = pb.Psi(
            0.0
        )

        Q_values[mask_near] = (
            0.5
            * np.sinh(psi_center)
            * R_array[mask_near]
        )

    # Regular expression away from the axis
    if np.any(mask_regular):

        R_regular = (
            R_array[mask_regular]
        )

        Q_values[mask_regular] = (
            hr.HR(R_regular)
            / (
                R_regular
                * mr.MR(R_regular)
            )
        )

    if np.ndim(R) == 0:
        return float(Q_values)

    return Q_values


# ============================================================
# AUXILIARY ODE
# ============================================================

def f1_ode_system(R, y):
    """
    Auxiliary ODE for F1(R):

        F1'(R) = Lambda delta^2 Q(R).
    """

    return [
        lam.Lambda
        * par.delta**2
        * QF1(R)
    ]


# ============================================================
# CACHED ODE SOLUTION
# ============================================================

_f1_solution = None


def get_f1_solution():
    """
    Return the cached auxiliary-ODE solution for F1(R).

    Since

        F1(1) = 0,

    solve_ivp integrates backward from R = 1 to R = 0.

    The ODE is solved only on the first call.
    """

    global _f1_solution

    if _f1_solution is None:

        _f1_solution = solve_ivp(
            f1_ode_system,
            (1.0, 0.0),
            [0.0],
            method="DOP853",
            rtol=RTOL_F1,
            atol=ATOL_F1,
            dense_output=True
        )

        if not _f1_solution.success:
            raise RuntimeError(
                "F1(R) auxiliary ODE solver failed: "
                + _f1_solution.message
            )

    return _f1_solution


# ============================================================
# PUBLIC FUNCTION: F1(R)
# ============================================================

def F1(R):
    """
    Dimensionless auxiliary function F1(R).

    Parameters
    ----------
    R : float or array_like
        Dimensionless radial coordinate.

        Required domain:
            0 <= R <= 1

    Returns
    -------
    float or ndarray
        F1(R).
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

    solution = get_f1_solution()

    F1_values = solution.sol(
        R_array
    )[0]

    # Enforce exact boundary condition:
    #
    #   F1(1) = 0.
    #
    F1_values = np.where(
        R_array == 1.0,
        0.0,
        F1_values
    )

    if np.ndim(R) == 0:
        return float(F1_values)

    return F1_values


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
        "F1(R) operational check"
    )

    print(
        "-----------------------"
    )

    print(
        f"F1(0)   = {F1(0.0):.15e}"
    )

    print(
        f"F1(0.5) = {F1(0.5):.15e}"
    )

    print(
        f"F1(1)   = {F1(1.0):.15e}"
    )

    print("")
    print(
        f"Q(0)    = {QF1(0.0):.15e}"
    )
