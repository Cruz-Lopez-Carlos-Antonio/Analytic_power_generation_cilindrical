import numpy as np
import matplotlib.pyplot as plt

from scipy.integrate import quad, solve_ivp

import Parameters as par
import Lambda_parameter as lam
import PoissonBoltzman as pb
import HR as hr
import MR as mr


# ============================================================
# F1(R) — VALIDATION
# ============================================================
#
# We compare
#
#   F1(R) =
#       -Lambda delta^2 Integral_R^1 Q(s) ds
#
# with the equivalent auxiliary ODE
#
#   F1'(R) = Lambda delta^2 Q(R)
#   F1(1)  = 0
#
# where
#
#   Q(R) = H(R)/(R M(R)).
#
# Near R = 0:
#
#   Q(R) =
#       [sinh(Psi(0))/2] R + O(R^3).
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

def Q(R):
    """
    Auxiliary function

        Q(R) = H(R)/(R M(R)).

    The apparent singularity at R = 0 is removable.

    Near the axis,

        Q(R) ~ 0.5*sinh(Psi(0))*R.
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

    # Exact removable value
    Q_values[mask_axis] = 0.0

    # Local asymptotic expansion
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
# METHOD 1: DIRECT QUADRATURE
# ============================================================

def F1Direct(
    R,
    epsabs=1.0e-12,
    epsrel=1.0e-12,
    limit=300
):
    """
    Direct quadrature:

        F1(R) =
            -Lambda delta^2
             Integral_R^1 Q(s) ds.
    """

    R = float(R)

    if not 0.0 <= R <= 1.0:
        raise ValueError(
            "R must satisfy 0 <= R <= 1."
        )

    if R == 1.0:
        return 0.0

    quad_kwargs = dict(
        epsabs=epsabs,
        epsrel=epsrel,
        limit=limit
    )

    # Tell QUADPACK about the transition point
    # only if it lies inside the integration interval.
    if R < EPS_Q:
        quad_kwargs["points"] = [
            EPS_Q
        ]

    integral_value, _ = quad(
        Q,
        R,
        1.0,
        **quad_kwargs
    )

    return (
        -lam.Lambda
        * par.delta**2
        * integral_value
    )


# ============================================================
# METHOD 2: AUXILIARY ODE
# ============================================================

def f1_ode_system(R, y):
    """
    Auxiliary ODE:

        F1'(R) = Lambda delta^2 Q(R).
    """

    return [
        lam.Lambda
        * par.delta**2
        * Q(R)
    ]


_f1_solution = None


def get_f1_solution():
    """
    Return the cached F1 auxiliary-ODE solution.

    Because the boundary condition is

        F1(1) = 0,

    solve_ivp integrates naturally backward from R=1 to R=0.
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
                "F1 auxiliary ODE solver failed: "
                + _f1_solution.message
            )

    return _f1_solution


def F1ODE(R):
    """
    F1(R) obtained from the auxiliary ODE.
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

    values = solution.sol(
        R_array
    )[0]

    # Exact boundary condition
    values = np.where(
        R_array == 1.0,
        0.0,
        values
    )

    if np.ndim(R) == 0:
        return float(values)

    return values


# ============================================================
# COMPARISON
# ============================================================

def validate_F1():

    R_test = np.linspace(
        0.0,
        1.0,
        20
    )

    direct_values = np.array([
        F1Direct(R)
        for R in R_test
    ])

    ode_values = F1ODE(
        R_test
    )

    abs_error = np.abs(
        direct_values
        - ode_values
    )

    ape = np.full_like(
        abs_error,
        np.nan
    )

    mask = (
        direct_values != 0.0
    )

    ape[mask] = (
        100.0
        * abs_error[mask]
        / np.abs(
            direct_values[mask]
        )
    )

    print("")
    print(
        "F1(R): direct quadrature vs auxiliary ODE"
    )

    print(
        "--------------------------------------------------------------"
    )

    print(
        "       R          F1Direct          F1ODE        Abs. error"
    )

    print(
        "--------------------------------------------------------------"
    )

    for R, direct, ode, error in zip(
        R_test,
        direct_values,
        ode_values,
        abs_error
    ):

        print(
            f"{R:8.6f}   "
            f"{direct: .15e}   "
            f"{ode: .15e}   "
            f"{error:.3e}"
        )

    print("")
    print(
        "Endpoint checks"
    )

    print(
        "---------------"
    )

    print(
        f"F1Direct(0) = "
        f"{F1Direct(0.0):.17e}"
    )

    print(
        f"F1ODE(0)    = "
        f"{F1ODE(0.0):.17e}"
    )

    print(
        f"F1Direct(1) = "
        f"{F1Direct(1.0):.17e}"
    )

    print(
        f"F1ODE(1)    = "
        f"{F1ODE(1.0):.17e}"
    )

    print("")
    print(
        "Axis regularity"
    )

    print(
        "---------------"
    )

    print(
        f"Q(0) = "
        f"{Q(0.0):.17e}"
    )

    print(
        "Expected F1'(0) = "
        f"{lam.Lambda * par.delta**2 * Q(0.0):.17e}"
    )

    print("")
    print(
        "Global comparison"
    )

    print(
        "-----------------"
    )

    print(
        "Maximum absolute error = "
        f"{np.max(abs_error):.15e}"
    )

    print(
        "Maximum APE (%)        = "
        f"{np.nanmax(ape):.15e}"
    )

    return (
        R_test,
        direct_values,
        ode_values,
        abs_error,
        ape
    )


# ============================================================
# PLOT
# ============================================================

def plot_F1():

    R_plot = np.linspace(
        0.0,
        1.0,
        500
    )

    ode_values = F1ODE(
        R_plot
    )

    direct_points = np.linspace(
        0.0,
        1.0,
        40
    )

    direct_values = np.array([
        F1Direct(R)
        for R in direct_points
    ])

    plt.figure()

    plt.plot(
        R_plot,
        ode_values,
        label="F1ODE"
    )

    plt.plot(
        direct_points,
        direct_values,
        "o",
        markersize=3,
        label="F1Direct"
    )

    plt.xlabel(
        r"$R$"
    )

    plt.ylabel(
        r"$F_1(R)$"
    )

    plt.title(
        r"$F_1(R)$: direct quadrature vs auxiliary ODE"
    )

    plt.legend()

    plt.tight_layout()

    plt.show()


# ============================================================
# DIRECT EXECUTION
# ============================================================

if __name__ == "__main__":

    validate_F1()

    plot_F1()
