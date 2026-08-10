# ============================================================
# Lambda_parameter.py
# ============================================================
#
# Computation of the dimensionless parameter
#
#                 2 epsilon^2 kappa^2 zeta_T^2
#       Lambda = -------------------------------
#                    mu sigma_inf F_cc
#
# The physical parameters are imported from Parameters.py,
# while F_cc is imported from the previously validated
# F_cc.py module.
#
# ============================================================


import Parameters as par
import F_cc as fcc


# ============================================================
# COMPUTE Lambda
# ============================================================

Lambda = (
    2.0
    * par.epsilon**2
    * par.kappa**2
    * par.zeta_T**2
    / (
        par.mu
        * par.sigma_inf
        * fcc.Fcc
    )
)


# ============================================================
# DIRECT EXECUTION
# ============================================================
#
# This block is executed only when this file is run directly.
#
# If another module uses
#
#       import Lambda_parameter as lam
#
# nothing is printed automatically.
#
# ============================================================

if __name__ == "__main__":

    print("")
    print("Lambda")
    print("------")

    print(
        f"Lambda = {Lambda:.17g}"
    )
