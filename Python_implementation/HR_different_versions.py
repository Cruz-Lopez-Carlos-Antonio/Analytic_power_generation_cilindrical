import numpy as np

from scipy.integrate import quad
from scipy.integrate import solve_ivp
from scipy.integrate import cumulative_trapezoid
from scipy.interpolate import PchipInterpolator

import PoissonBoltzman as pb


# ============================================================
# H(R) — DIFFERENT NUMERICAL VERSIONS
# ============================================================
#
# We define
#
#   H(R) = integral_0^R tau*sinh(Psi(tau)) dtau
#
# and compare three numerical approaches:
#
#   1. Direct quadrature with scipy.integrate.quad
#   2. Auxiliary ODE:
#
#          H'(R) = R*sinh(Psi(R)),
#          H(0)  = 0
#
#   3. Cumulative trapezoidal quadrature + interpolation
#
# ============================================================


# ============================================================
# COMMON INTEGRAND
# ============================================================

def h_integrand(R):
    """
    Integrand of H(R):

        R*sinh(Psi(R)).
    """

    return R * np.sinh(pb.Psi(R))


# ============================================================
# METHOD 1 — DIRECT QUADRATURE
# ============================================================

def H_direct(R):
    """
    Direct evaluation of

        H(R) = integral_0^R tau*sinh(Psi(tau)) dtau

    using scipy.integrate.quad.
    """

    R = float(R)

    if R < 0.0 or R > 1.0:
        raise ValueError(
            "R must satisfy 0 <= R <= 1."
        )

    if R == 0.0:
        return 0.0

    value, error = quad(
        h_integrand,
        0.0,
        R,
        epsabs=1.0e-13,
        epsrel=1.0e-11,
        limit=300
    )

    return value


# ============================================================
# METHOD 2 — AUXILIARY ODE
# ============================================================

def h_ode_system(R, y):
    """
    Auxiliary ODE:

        H'(R) = R*sinh(Psi(R)).
    """

    return [
        R * np.sinh(pb.Psi(R))
    ]


_h_ode_solution = None


def get_h_ode_solution():
    """
    Solve the auxiliary ODE once and cache
    the dense-output solution.
    """

    global _h_ode_solution

    if _h_ode_solution is None:

        _h_ode_solution = solve_ivp(
            h_ode_system,
            (0.0, 1.0),
            [0.0],
            method="DOP853",
            rtol=1.0e-11,
            atol=1.0e-13,
            dense_output=True
        )

        if not _h_ode_solution.success:
            raise RuntimeError(
                "Auxiliary H(R) ODE solver failed: "
                + _h_ode_solution.message
            )

    return _h_ode_solution


def H_ode(R):
    """
    Evaluate H(R) from the auxiliary ODE.
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

    solution = get_h_ode_solution()

    values = solution.sol(
        R_array
    )[0]

    values = np.where(
        R_array == 0.0,
        0.0,
        values
    )

    if np.ndim(R) == 0:
        return float(values)

    return values


# ============================================================
# METHOD 3 — CUMULATIVE TRAPEZOID + INTERPOLATION
# ============================================================

_h_trap_interpolator = None


def get_h_trap_interpolator(
    n_grid=5001
):
    """
    Construct a cumulative trapezoidal approximation
    of H(R) on a uniform grid and interpolate it.
    """

    global _h_trap_interpolator

    if _h_trap_interpolator is None:

        R_grid = np.linspace(
            0.0,
            1.0,
            n_grid
        )

        integrand_grid = (
            R_grid
            * np.sinh(
                pb.Psi(R_grid)
            )
        )

        H_grid = cumulative_trapezoid(
            integrand_grid,
            R_grid,
            initial=0.0
        )

        _h_trap_interpolator = PchipInterpolator(
            R_grid,
            H_grid
        )

    return _h_trap_interpolator


def H_trap(R):
    """
    Evaluate H(R) using cumulative trapezoidal
    quadrature and PCHIP interpolation.
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

    interpolator = get_h_trap_interpolator()

    values = interpolator(
        R_array
    )

    values = np.where(
        R_array == 0.0,
        0.0,
        values
    )

    if np.ndim(R) == 0:
        return float(values)

    return values


# ============================================================
# 10-POINT CROSS-VALIDATION
# ============================================================

def compare_methods():
    """
    Compare the three methods at 10 uniformly
    spaced radial positions.

    H_direct is used as the provisional reference.
    """

    R_sample = np.linspace(
        0.0,
        1.0,
        10
    )

    print(
        "\nH(R) comparison"
    )

    print(
        "-----------------------------------------------------------------------------------------------"
    )

    print(
        "       R            H direct            H ODE"
        "              H trap          APE ODE (%)       APE trap (%)"
    )

    print(
        "-----------------------------------------------------------------------------------------------"
    )

    for R in R_sample:

        hd = H_direct(R)
        ho = H_ode(R)
        ht = H_trap(R)

        if R == 0.0:

            ape_ode = np.nan
            ape_trap = np.nan

        else:

            ape_ode = (
                100.0
                * abs(
                    (ho - hd) / hd
                )
            )

            ape_trap = (
                100.0
                * abs(
                    (ht - hd) / hd
                )
            )

        print(
            f"{R: .10f}   "
            f"{hd: .15e}   "
            f"{ho: .15e}   "
            f"{ht: .15e}   "
            f"{ape_ode: .8e}   "
            f"{ape_trap: .8e}"
        )


# ============================================================
# ABSOLUTE-ERROR CHECK
# ============================================================

def compare_absolute_errors():
    """
    Compare absolute errors against H_direct.
    """

    R_sample = np.linspace(
        0.0,
        1.0,
        10
    )

    print(
        "\nAbsolute errors"
    )

    print(
        "--------------------------------------------------------------"
    )

    print(
        "       R          |ODE-direct|          |trap-direct|"
    )

    print(
        "--------------------------------------------------------------"
    )

    for R in R_sample:

        hd = H_direct(R)
        ho = H_ode(R)
        ht = H_trap(R)

        print(
            f"{R: .10f}   "
            f"{abs(ho - hd): .15e}   "
            f"{abs(ht - hd): .15e}"
        )


# ============================================================
# ASYMPTOTIC AXIS CHECK
# ============================================================

def validate_axis_behavior():
    """
    Check

        H(R)/R^2 -> sinh(Psi(0))/2

    as R -> 0.
    """

    R_sample = np.array([
        1.0e-3,
        2.0e-3,
        5.0e-3,
        1.0e-2
    ])

    theoretical_limit = (
        np.sinh(
            pb.Psi(0.0)
        )
        / 2.0
    )

    print(
        "\nAxis asymptotic check"
    )

    print(
        "-----------------------------------------------------------------------"
    )

    print(
        "       R              H_ode(R)/R^2"
        "          sinh(Psi(0))/2        abs. difference"
    )

    print(
        "-----------------------------------------------------------------------"
    )

    for R in R_sample:

        ratio = (
            H_ode(R)
            / R**2
        )

        difference = abs(
            ratio
            - theoretical_limit
        )

        print(
            f"{R: .6e}   "
            f"{ratio: .15e}   "
            f"{theoretical_limit: .15e}   "
            f"{difference: .15e}"
        )


# ============================================================
# DIRECT EXECUTION
# ============================================================
#
# Nothing below is executed when this module is imported.
#
# ============================================================

if __name__ == "__main__":

    compare_methods()

    compare_absolute_errors()

    validate_axis_behavior()
