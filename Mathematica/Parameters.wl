(* ========================================================= *)
(* PARAMETERS                                                *)
(* ========================================================= *)
(*
   Physical constants, physical parameters, derived
   quantities and dimensionless parameters.

   Numerical values are stored with arbitrary precision
   so that downstream solvers such as NDSolveValue and
   FindRoot can work consistently with WorkingPrecision
   greater than MachinePrecision.

   IMPORTANT:
   Increased computational precision does NOT imply
   increased physical or experimental accuracy.
*)
(* ========================================================= *)


(* ========================================================= *)
(* 1. WORKING PRECISION FOR PARAMETERS                       *)
(* ========================================================= *)

wpPar = 50;

(* ========================================================= *)
(* 2. FUNDAMENTAL PHYSICAL CONSTANTS                         *)
(* ========================================================= *)

(* Boltzmann constant [J/K] *)
kB = N[1380649/10^29, wpPar];

(* Elementary charge [C] *)
e = N[1602176634/10^28, wpPar];

(* Avogadro constant [1/mol] *)
NA = N[602214076/10^-15, wpPar];


(* ========================================================= *)
(* 3. PHYSICAL PARAMETERS OF THE PROBLEM                     *)
(* ========================================================= *)

(* Temperature [K] *)
T = N[298, wpPar];

(* Fluid density [kg/m^3] *)
rho = N[1000, wpPar];

(* Dynamic viscosity [Pa s] *)
mu = N[891/10^6, wpPar];

(* Microchannel length [m] *)
l = N[1/10^3, wpPar];

(* Microchannel radius [m] *)
a = N[1/10^8, wpPar];

(* Characteristic water velocity [m/s] *)
Jw = N[868327/10^12, wpPar];

(* Dielectric permittivity [C V^-1 m^-1] = [F/m] *)
epsilon = N[6954/10^13, wpPar];

(* Ionic valence *)
z = N[1, wpPar];

(* Bulk molar concentration [mol/L] *)
Cbulk = N[1/2, wpPar];

(* Material parameter appearing in the viscosity model
   [m^2/V^2] *)
f = N[23/10^17, wpPar];

(* Ionic diffusion coefficient [m^2/s] *)
DiffCoeff = N[1312/10^12, wpPar];

(* Zeta potential at the wall [V] *)
zeta = N[-3165/10^5, wpPar];


(* ========================================================= *)
(* 4. DERIVED PHYSICAL QUANTITIES                            *)
(* ========================================================= *)

(* Number density [m^-3]
   nInf = 1000 Cbulk NA
*)
nInf =1000*Cbulk*NA;

(* Inverse Debye length kappa [m^-1]
   kappa =Sqrt[2 z^2 e^2 nInf /
      (epsilon kB T)]
*)

kappa =Sqrt[(2*z^2*e^2*nInf)/(epsilon*kB*T)];

(* Debye length [m]
   lambdaD = 1/kappa
*)

lambdaD =1/kappa;
(* Thermal potential [V]
   zetaT = kB T / (z e)
*)
zetaT =kB*T/(z*e);

(* Bulk osmotic pressure of the draw solution [Pa]

   Van't Hoff equation for NaCl:

   PosDb = 2 kB T nInf
*)

PosDb =2*kB*T*nInf;
(* Bulk electrical conductivity [S/m]
   sigmaInf =2 e^2 z^2 DiffCoeff nInf / (kB T)
*)

sigmaInf =2*e^2*z^2*DiffCoeff*nInf/(kB*T);

(* Dimensionless surface potential
   PsiS = zeta / zetaT
*)

PsiS =zeta/zetaT;

(* ========================================================= *)
(* 5. DIMENSIONLESS PARAMETERS                               *)
(* ========================================================= *)

(* Reynolds number
   Re = rho Jw l / mu
*)
Reynolds =rho*Jw*l/mu;

(* Geometric aspect ratio
  xi = a/l*)

xi =a/l;

(* Dimensionless Debye parameter
   delta = a kappa
         = a/lambdaD
*)

delta =a*kappa;

(* Dimensionless parameter associated with the
   electric-field-dependent viscosity model

   omega = f (zetaT/a)^2
*)

omega =f*(zetaT/a)^2;


(* Ratio of osmotic pressure to the characteristic
   pressure associated with the flow

   alpha = PosDb l / (mu Jw)
*)

alpha =
  PosDb*l/(mu*Jw);

