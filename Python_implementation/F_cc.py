import numpy as np

from scipy.integrate import quad

import PoissonBoltzman as pb


# ============================================================
# F_cc
# ============================================================
#
# We compute
#
#       F_cc = 2 Integral_0^1 R cosh(Psi(R)) dR
#
# using the previously validated Poisson-Boltzmann
# solution Psi(R).
#
# ============================================================


# ============================================================
# INTEGRAND
# ============================================================

def fcc_integrand(R):
    """
    Integrand appearing in

        F_cc = 2 Integral_0^1 R cosh(Psi(R)) dR.
    """

    return (
        2.0
        * R
        * np.cosh(pb.Psi(R))
    )


# ============================================================
# COMPUTE F_cc
# ============================================================

def compute_Fcc(
    epsabs=1.0e-13,
    epsrel=1.0e-13,
    limit=200
):
    """
    Compute F_cc by adaptive numerical quadrature.

    Parameters
    ----------
    epsabs : float
        Absolute error tolerance passed to scipy.integrate.quad.

    epsrel : float
        Relative error tolerance passed to scipy.integrate.quad.

    limit : int
        Maximum number of subintervals used by QUADPACK.

    Returns
    -------
    Fcc : float
        Numerical value of F_cc.

    error_estimate : float
        Absolute error estimate returned by quad.
    """

    Fcc, error_estimate = quad(
        fcc_integrand,
        0.0,
        1.0,
        epsabs=epsabs,
        epsrel=epsrel,
        limit=limit
    )

    return Fcc, error_estimate


# ============================================================
# MODULE-LEVEL VALUE
# ============================================================

Fcc, Fcc_error = compute_Fcc()


