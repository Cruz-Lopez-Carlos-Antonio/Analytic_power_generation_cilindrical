import numpy as np
import matplotlib.pyplot as plt

from scipy.integrate import quad, solve_ivp

import Parameters as par
import MR as mr


# ============================================================
# F0(R) — VALIDATION
# ============================================================
#
# We compare
#
#   F0(R) =
#       -Integral_R^1 [PiD s/(2 M(s))] ds
#
# with the equivalent auxiliary ODE
#
#   F0'(R) = PiD R/(2 M(R))
#   F0(1)  = 0
#
# where
#
#   PiD = alpha xi^2.
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
# METHOD 1: DIRECT QUADRATURE
# ============================================================

def f0_integrand(s):
    """
    Integrand appearing in

        F0(R) =
            -Integral_R^1 PiD*s/(2 M(s)) ds.
    """

    return (
        PiD
        * s
        / (
            2.0
            * mr.MR(s)
        )
    )


def F0Direct(
    R,
    epsabs=1.0e-11,
    epsrel=1.0e-11,
    limit=300
):
    """
    Compute F0(R) by direct adaptive quadrature.
    """

    R = float(R)

    if not 0.0 <= R <= 1.0:
        raise ValueError(
            "R must satisfy 0 <= R <= 1."
        )

    if R == 1.0:
        return 0.0

    integral_value, _ = quad(
        f0_integrand,
        R,
        1.0,
        epsabs=epsabs,
        epsrel=epsrel,
        limit=limit
    )

    return -integral_value


# ============================================================
# METHOD 2: AUXILIARY ODE
# ============================================================

def f0_ode_system(R, y):
    """
    Auxiliary ODE:

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


_f0_solution = None


def get_f0_solution():
    """
    Return the cached auxiliary-ODE solution.

    Since

        F0(1) = 0,

    solve_ivp integrates backward from R = 1 to R = 0.
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
                "F0 auxiliary ODE solver failed: "
                + _f0_solution.message
            )

    return _f0_solution


def F0ODE(R):
    """
    F0(R) obtained from the auxiliary ODE.
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

    values = solution.sol(
        R_array
    )[0]

    # Enforce exact boundary condition
    #
    #   F0(1) = 0.
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
# THEORETICAL DERIVATIVE
# ============================================================

def F0Prime(R):
    """
    Exact right-hand side of the auxiliary ODE:

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
# NUMERICAL DERIVATIVE OF ODE SOLUTION
# ============================================================

def F0PrimeNumerical(R, h=1.0e-6):
    """
    Numerical derivative of the dense ODE solution.

    A centered difference is used away from the endpoints.
    """

    R = float(R)

    if not 0.0 < R < 1.0:
        raise ValueError(
            "Numerical derivative requires 0 < R < 1."
        )

    h_eff = min(
        h,
        0.5 * R,
        0.5 * (1.0 - R)
    )

    return (
        F0ODE(R + h_eff)
        - F0ODE(R - h_eff)
    ) / (
        2.0 * h_eff
    )


# ============================================================
# VALIDATION
# ============================================================

def validate_F0():

    R_test = np.linspace(
        0.0,
        1.0,
        20
    )

    direct_values = np.array([
        F0Direct(R)
        for R in R_test
    ])

    ode_values = F0ODE(
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
        "F0(R): direct quadrature vs auxiliary ODE"
    )

    print(
        "---------------------------------------------------------------------"
    )

    print(
        "       R          F0Direct          F0ODE        Abs. error      APE (%)"
    )

    print(
        "---------------------------------------------------------------------"
    )

    for R, direct, ode, error, ape_i in zip(
        R_test,
        direct_values,
        ode_values,
        abs_error,
        ape
    ):

        ape_text = (
            "       ---"
            if np.isnan(ape_i)
            else f"{ape_i:.3e}"
        )

        print(
            f"{R:8.6f}   "
            f"{direct: .15e}   "
            f"{ode: .15e}   "
            f"{error:.3e}   "
            f"{ape_text}"
        )

    # --------------------------------------------------------
    # Endpoint checks
    # --------------------------------------------------------

    print("")
    print(
        "Endpoint checks"
    )

    print(
        "---------------"
    )

    print(
        f"PiD         = {PiD:.17e}"
    )

    print(
        f"F0Direct(0) = {F0Direct(0.0):.17e}"
    )

    print(
        f"F0ODE(0)    = {F0ODE(0.0):.17e}"
    )

    print(
        f"F0Direct(1) = {F0Direct(1.0):.17e}"
    )

    print(
        f"F0ODE(1)    = {F0ODE(1.0):.17e}"
    )

    # --------------------------------------------------------
    # Axis regularity
    # --------------------------------------------------------

    print("")
    print(
        "Axis regularity"
    )

    print(
        "---------------"
    )

    print(
        f"M(0)            = {mr.MR(0.0):.17e}"
    )

    print(
        f"Expected F0'(0) = {F0Prime(0.0):.17e}"
    )

    # --------------------------------------------------------
    # Differential residual
    # --------------------------------------------------------

    residual_points = np.linspace(
        0.05,
        0.95,
        19
    )

    residuals = np.array([
        F0PrimeNumerical(R)
        - F0Prime(R)
        for R in residual_points
    ])

    print("")
    print(
        "Differential residual"
    )

    print(
        "---------------------"
    )

    print(
        "Maximum |F0'_num(R) - PiD R/(2 M(R))| = "
        f"{np.max(np.abs(residuals)):.15e}"
    )

    # --------------------------------------------------------
    # Global comparison
    # --------------------------------------------------------

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

    # --------------------------------------------------------
    # Sign check
    # --------------------------------------------------------

    sign_check = np.all(
        ode_values[:-1] < 0.0
    )

    print("")
    print(
        "Sign check"
    )

    print(
        "----------"
    )

    print(
        "F0(R) < 0 for 0 <= R < 1 : "
        f"{sign_check}"
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

def plot_F0():

    R_plot = np.linspace(
        0.0,
        1.0,
        500
    )

    ode_values = F0ODE(
        R_plot
    )

    direct_points = np.linspace(
        0.0,
        1.0,
        40
    )

    direct_values = np.array([
        F0Direct(R)
        for R in direct_points
    ])

    plt.figure()

    plt.plot(
        R_plot,
        ode_values,
        label="F0ODE"
    )

    plt.plot(
        direct_points,
        direct_values,
        "o",
        markersize=3,
        label="F0Direct"
    )

    plt.xlabel(
        r"$R$"
    )

    plt.ylabel(
        r"$F_0(R)$"
    )

    plt.title(
        r"$F_0(R)$: direct quadrature vs auxiliary ODE"
    )

    plt.legend()

    plt.tight_layout()

    plt.show()


# ============================================================
# DIRECT EXECUTION
# ============================================================

if __name__ == "__main__":

    validate_F0()

    plot_F0()
