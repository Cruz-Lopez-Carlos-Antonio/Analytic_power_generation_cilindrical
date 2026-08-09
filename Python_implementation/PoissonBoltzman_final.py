import numpy as np
import matplotlib.pyplot as plt

from scipy.integrate import solve_bvp
from scipy.special import iv

import Parameters as par


# ============================================================
# POISSON-BOLTZMANN EQUATION
# ============================================================
#
# We solve
#
#   Psi''(R) + (1/R) Psi'(R) = delta^2 sinh(Psi(R))
#
# for
#
#   0 <= R <= 1
#
# with boundary conditions
#
#   Psi'(0) = 0
#   Psi(1)  = Psi_s
#
#
# To avoid the apparent singularity at R = 0, define
#
#   p(R) = R Psi'(R)
#
# so that
#
#   Psi' = p/R
#   p'   = delta^2 R sinh(Psi)
#
# and the boundary conditions become
#
#   p(0)   = 0
#   Psi(1) = Psi_s
#
# ============================================================


def pb_system(R, y):
    """
    First-order form of the cylindrical Poisson-Boltzmann equation.

    y[0] = Psi(R)
    y[1] = p(R) = R*Psi'(R)
    """

    Psi_value = y[0]
    p = y[1]

    # Safe evaluation of p/R.
    #
    # At R = 0, regularity implies
    #
    #   p(R) = O(R^2),
    #
    # hence
    #
    #   Psi'(0) = 0.
    #
    dPsi_dR = np.divide(
        p,
        R,
        out=np.zeros_like(p),
        where=(R != 0.0)
    )

    dp_dR = (
        par.delta**2
        * R
        * np.sinh(Psi_value)
    )

    return np.vstack(
        (dPsi_dR, dp_dR)
    )


# ============================================================
# BOUNDARY CONDITIONS
# ============================================================

def pb_boundary_conditions(ya, yb):
    """
    Boundary conditions:

        p(0)   = 0
        Psi(1) = Psi_s
    """

    return np.array([
        ya[1],
        yb[0] - par.Psi_s
    ])


# ============================================================
# INITIAL GUESS
# ============================================================

def initial_guess(R):
    """
    Initial guess based on the linearized
    Poisson-Boltzmann equation.

    For small Psi,

        sinh(Psi) ~ Psi,

    giving the cylindrical Debye-Huckel solution

        Psi(R)
        ~ Psi_s I0(delta R) / I0(delta).

    This approximation is used ONLY as an
    initial guess for the nonlinear BVP solver.
    """

    denominator = iv(
        0,
        par.delta
    )

    Psi_guess = (
        par.Psi_s
        * iv(0, par.delta * R)
        / denominator
    )

    # Since
    #
    #   p(R) = R Psi'(R),
    #
    # the corresponding initial guess is
    #
    p_guess = (
        par.Psi_s
        * par.delta
        * R
        * iv(1, par.delta * R)
        / denominator
    )

    return np.vstack(
        (Psi_guess, p_guess)
    )


# ============================================================
# BVP SOLVER
# ============================================================

def solve_poisson_boltzmann(
    n_mesh=400,
    tol=1.0e-8,
    max_nodes=20000
):
    """
    Solve the nonlinear Poisson-Boltzmann
    boundary-value problem.
    """

    R = np.linspace(
        0.0,
        1.0,
        n_mesh
    )

    y_guess = initial_guess(R)

    solution = solve_bvp(
        pb_system,
        pb_boundary_conditions,
        R,
        y_guess,
        tol=tol,
        max_nodes=max_nodes,
        verbose=0
    )

    if not solution.success:
        raise RuntimeError(
            "Poisson-Boltzmann solver failed: "
            + solution.message
        )

    return solution


# ============================================================
# CACHED MODULE-LEVEL SOLUTION
# ============================================================
#
# The BVP is solved only when it is first needed.
# The resulting solution is then stored and reused.
#
# Thus:
#
#   Psi(0.5)
#   Psi(0.75)
#   PsiPrime(0.5)
#
# do NOT solve the BVP three times.
#
# ============================================================

_pb_solution = None


def get_solution():
    """
    Return the cached Poisson-Boltzmann solution.

    The BVP is solved only on the first call.
    """

    global _pb_solution

    if _pb_solution is None:
        _pb_solution = solve_poisson_boltzmann()

    return _pb_solution


# ============================================================
# PUBLIC FUNCTION: Psi(R)
# ============================================================

def Psi(R):
    """
    Dimensionless electrostatic potential Psi(R).

    Parameters
    ----------
    R : float or array_like
        Dimensionless radial coordinate.

        Required domain:
            0 <= R <= 1

    Returns
    -------
    float or ndarray
        Psi(R).
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

    solution = get_solution()

    Psi_values = solution.sol(
        R_array
    )[0]

    # Preserve scalar input -> scalar output
    if np.ndim(R) == 0:
        return float(Psi_values)

    return Psi_values


# ============================================================
# PUBLIC FUNCTION: p(R) = R Psi'(R)
# ============================================================

def pPB(R):
    """
    Auxiliary variable

        p(R) = R Psi'(R).
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

    solution = get_solution()

    p_values = solution.sol(
        R_array
    )[1]

    if np.ndim(R) == 0:
        return float(p_values)

    return p_values


# ============================================================
# PUBLIC FUNCTION: Psi'(R)
# ============================================================

def PsiPrime(R):
    """
    Radial derivative Psi'(R).

    Regularity at the symmetry axis gives

        Psi'(0) = 0.
    """

    R_array = np.asarray(
        R,
        dtype=float
    )

    p_values = np.asarray(
        pPB(R_array),
        dtype=float
    )

    derivative = np.divide(
        p_values,
        R_array,
        out=np.zeros_like(
            R_array,
            dtype=float
        ),
        where=(R_array != 0.0)
    )

    if np.ndim(R) == 0:
        return float(derivative)

    return derivative


# ============================================================
# DENSE EVALUATION
# ============================================================

def evaluate_solution(
    solution=None,
    n_points=1000
):
    """
    Evaluate Psi(R), Psi'(R), and p(R)
    on a dense uniform mesh.
    """

    if solution is None:
        solution = get_solution()

    R = np.linspace(
        0.0,
        1.0,
        n_points
    )

    values = solution.sol(R)

    Psi_values = values[0]
    p_values = values[1]

    dPsi_dR = np.divide(
        p_values,
        R,
        out=np.zeros_like(p_values),
        where=(R != 0.0)
    )

    return (
        R,
        Psi_values,
        dPsi_dR,
        p_values
    )


# ============================================================
# VALIDATION
# ============================================================

def validate_solution(solution=None):
    """
    Basic numerical checks for the
    Poisson-Boltzmann solution.
    """

    if solution is None:
        solution = get_solution()

    (
        R,
        Psi_values,
        dPsi_dR,
        p_values
    ) = evaluate_solution(solution)

    print(
        "\nPoisson-Boltzmann validation"
    )

    print(
        "----------------------------"
    )

    print(
        f"Solver success       = "
        f"{solution.success}"
    )

    print(
        f"Solver message       = "
        f"{solution.message}"
    )

    print(
        f"Psi(0)               = "
        f"{Psi_values[0]:.15e}"
    )

    print(
        f"Psi'(0)              = "
        f"{dPsi_dR[0]:.15e}"
    )

    print(
        f"Psi(1)               = "
        f"{Psi_values[-1]:.15e}"
    )

    print(
        f"Psi_s                = "
        f"{par.Psi_s:.15e}"
    )

    print(
        f"|Psi(1)-Psi_s|       = "
        f"{abs(Psi_values[-1] - par.Psi_s):.15e}"
    )

    print(
        f"p(0)                 = "
        f"{p_values[0]:.15e}"
    )

    if solution.rms_residuals.size > 0:

        print(
            f"Maximum RMS residual = "
            f"{np.max(solution.rms_residuals):.15e}"
        )


# ============================================================
# SELECTED VALUES FOR CROSS-VALIDATION
# ============================================================

def sample_solution(solution=None):
    """
    Evaluate Psi(R) at selected radial positions
    for cross-validation.
    """

    if solution is None:
        solution = get_solution()

    R_sample = np.array([
        0.0,
        0.25,
        0.50,
        0.75,
        1.0
    ])

    Psi_sample = solution.sol(
        R_sample
    )[0]

    print(
        "\nSelected Poisson-Boltzmann values"
    )

    print(
        "---------------------------------"
    )

    print(
        "       R                  Psi(R)"
    )

    print(
        "---------------------------------"
    )

    for R_i, Psi_i in zip(
        R_sample,
        Psi_sample
    ):
        print(
            f"{R_i:8.2f}     "
            f"{Psi_i: .15e}"
        )

    return (
        R_sample,
        Psi_sample
    )


# ============================================================
# PLOT
# ============================================================

def plot_solution(solution=None):
    """
    Plot the dimensionless EDL potential Psi(R).
    """

    if solution is None:
        solution = get_solution()

    (
        R,
        Psi_values,
        _,
        _
    ) = evaluate_solution(solution)

    plt.figure()

    plt.plot(
        R,
        Psi_values
    )

    plt.xlabel(
        r"$R$"
    )

    plt.ylabel(
        r"$\Psi(R)$"
    )

    plt.title(
        "Poisson-Boltzmann potential"
    )

    plt.grid(True)

    plt.tight_layout()

    plt.show()


# ============================================================
# DIRECT EXECUTION
# ============================================================
#
# Everything below this condition is executed only when
# this file itself is run directly.
#
# It is NOT executed when another script imports
#
#   PoissonBoltzman_final
#
# ============================================================

if __name__ == "__main__":

    solution = get_solution()

    validate_solution(
        solution
    )

    sample_solution(
        solution
    )

    plot_solution(
        solution
    )
