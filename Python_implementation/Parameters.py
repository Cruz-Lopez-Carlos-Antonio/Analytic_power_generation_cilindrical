import numpy as np


# ============================================================
# FUNDAMENTAL PHYSICAL CONSTANTS
# ============================================================

# Boltzmann constant [J/K]
k_B = 1.380649e-23

# Elementary charge [C]
e = 1.602176634e-19

# Avogadro constant [1/mol]
N_A = 6.02214076e23


# ============================================================
# PHYSICAL PARAMETERS OF THE PROBLEM
# ============================================================

# Temperature [K]
T = 298.0

# Fluid density [kg/m^3]
rho = 1000.0

# Dynamic viscosity [Pa s]
mu = 0.891e-3

# Microchannel length [m]
l = 1.0e-3

# Microchannel radius [m]
a = 1.0e-8

# Characteristic water velocity [m/s]
J_w = 8.68327e-7

# Dielectric permittivity of the electrolyte
# [C V^-1 m^-1] = [F/m]
epsilon = 6.954e-10

# Ionic valence
z = 1.0

# Bulk molar concentration [mol/L]
C_bulk = 0.5

# Material parameter appearing in the viscosity model [m^2/V^2]
f = 2.3e-16

# Ionic diffusion coefficient [m^2/s]
D = 1.312e-9

# Zeta potential at the wall [V]
zeta = -0.03165

# ============================================================
# DERIVED PHYSICAL QUANTITIES
# ============================================================

# Number density of one ionic species [m^-3]
#
# Conversion:
#   C_bulk [mol/L]
#       -> 1000*C_bulk [mol/m^3]
#       -> 1000*C_bulk*N_A [particles/m^3]
#
n_inf = 1000.0 * C_bulk * N_A


# Inverse Debye length, kappa [m^-1]
#
# Definition:
#   kappa^(-1) =
#       sqrt(epsilon*k_B*T /
#            (2*z^2*e^2*n_inf))
# therefore
#   kappa =
#       sqrt(2*z^2*e^2*n_inf /
#            (epsilon*k_B*T))
#
kappa = np.sqrt(2.0 * z**2 * e**2 * n_inf/ (epsilon * k_B * T))

# Debye length [m]
# lambda_D = 1/kappa

lambda_D = 1.0 / kappa

# Thermal potential [V]
#
#   zeta_T = k_B*T/(z*e)
#
zeta_T = k_B * T / (z * e)


# Bulk osmotic pressure of the draw solution [Pa]
#
# Van't Hoff equation for NaCl (i = 2):
#
#   P_osD_b = 2*k_B*T*n_inf
#
# Equivalent to:
#
#   P_osD_b = 2*R*T*C_bulk
#
# when C_bulk is expressed in mol/m^3.
#
P_osD_b = 2.0 * k_B * T * n_inf
# Bulk electrical conductivity [S/m]
#
#   sigma_inf = 2*e^2*z^2*D*n_inf/(k_B*T)
#
sigma_inf = (2.0 * e**2 * z**2 * D * n_inf/ (k_B * T))

# Dimensionless surface potential
#
#   Psi_s = zeta / zeta_T
#
Psi_s = zeta / zeta_T

# ============================================================
# DIMENSIONLESS PARAMETERS
# ============================================================

# Reynolds number
#
#   Re = rho*J_w*l/mu
#
# Ratio between inertial and viscous effects,
# using J_w as the characteristic velocity.
#
Re = rho * J_w * l / mu


# Geometric aspect ratio
#
#   xi = a/l
#
# Ratio between microchannel radius and length.
#
xi = a / l


# Dimensionless Debye parameter
#
#   delta = a*kappa = a/lambda_D
#
# Ratio between microchannel radius and Debye length.
#
delta = a * kappa

# Dimensionless parameter associated with the
# electric-field-dependent viscosity model
#
#   omega = f*(zeta_T/a)^2
#
omega = f * (zeta_T / a)**2

# Ratio of osmotic pressure to the characteristic
# pressure associated with the flow
#
#   alpha = P_osD_b*l/(mu*J_w)
#
alpha = P_osD_b * l / (mu * J_w)




# ============================================================
# OUTPUT
# ============================================================
if __name__ == "__main__":
    print("Derived physical quantities")
    print("---------------------------")
    print(f"n_inf     = {n_inf:.15e} m^-3")
    print(f"kappa     = {kappa:.15e} m^-1")
    print(f"lambda_D  = {lambda_D:.15e} m")
    print(f"zeta_T    = {zeta_T:.15e} V")
    print(f"P_osD_b   = {P_osD_b:.15e} Pa")
    print(f"sigma_inf = {sigma_inf:.15e} S/m")

    print("\nDimensionless parameters")
    print("------------------------")
    print(f"Re        = {Re:.15e}")
    print(f"xi        = {xi:.15e}")
    print(f"delta     = {delta:.15e}")
    print(f"omega     = {omega:.15e}")
    print(f"alpha     = {alpha:.15e}")
