(* ============================================================ *)
(* Phi_validate.wl                                              *)
(* Verification of the dimensionless streaming potential        *)
(* ============================================================ *)

ClearAll["Global`*"];

Module[{localDir = DirectoryName[$InputFileName]},

  (* Previous validated modules *)
  Get[FileNameJoin[{localDir, "G.wl"}]];
  Get[FileNameJoin[{localDir, "Omega_parameter.wl"}]];
  Get[FileNameJoin[{localDir, "Lambda_parameter.wl"}]];
];

(* ============================================================ *)
(* NUMERICAL SETTINGS                                           *)
(* ============================================================ *)

wpPhi = 30;
agPhi = 12;
pgPhi = 12;

(* ============================================================ *)
(* STREAMING-POTENTIAL GRADIENT                                 *)
(* ============================================================ *)
(*                                                              *)
(* From Eq. (56):                                               *)
(*                                                              *)
(*   dPhi/dZ = -Lambda delta^2 G Omega                          *)
(*                                                              *)
(* ============================================================ *)

ClearAll[PhiGradient];

PhiGradient =
  -Lambda*delta^2*G*Omega;

(* ============================================================ *)
(* STREAMING POTENTIAL                                          *)
(* ============================================================ *)
(*                                                              *)
(* With Phi(0) = 0:                                             *)
(*                                                              *)
(*   Phi(Z) = PhiGradient Z                                     *)
(*                                                              *)
(* ============================================================ *)

ClearAll[Phi];

Phi[Z_?NumericQ] := PhiGradient*Z;

(* ============================================================ *)
(* AVERAGE STREAMING POTENTIAL                                  *)
(* ============================================================ *)
(*                                                              *)
(* Eq. (58):                                                    *)
(*                                                              *)
(*   PhiAv = Integral_0^1 Phi(Z) dZ                             *)
(*         = PhiGradient/2                                      *)
(*                                                              *)
(* ============================================================ *)

ClearAll[PhiAvAnalytical, PhiAvNumerical];

PhiAvAnalytical = PhiGradient/2;

PhiAvNumerical = NIntegrate[
   Phi[Z],
   {Z, 0, 1},
   WorkingPrecision -> wpPhi,
   AccuracyGoal -> agPhi,
   PrecisionGoal -> pgPhi,
   MaxRecursion -> 30,
   Method -> "GlobalAdaptive"
];

(* ============================================================ *)
(* VERIFICATION                                                 *)
(* ============================================================ *)

Print["================================================"];
Print["Streaming potential verification"];
Print["================================================"];

Print["delta          = ", N[delta, 18]];
Print["Lambda         = ", N[Lambda, 18]];
Print["Omega          = ", N[Omega, 18]];
Print["G              = ", N[G, 18]];

Print[""];
Print["dPhi/dZ        = ", N[PhiGradient, 18]];

Print[""];
Print["Phi(0)         = ", N[Phi[0], 18]];
Print["Phi(0.25)      = ", N[Phi[0.25], 18]];
Print["Phi(0.50)      = ", N[Phi[0.50], 18]];
Print["Phi(0.75)      = ", N[Phi[0.75], 18]];
Print["Phi(1)         = ", N[Phi[1], 18]];

Print[""];
Print["PhiAv analytic = ", N[PhiAvAnalytical, 18]];
Print["PhiAv numeric  = ", N[PhiAvNumerical, 18]];

Print[
  "Absolute residual = ",
  N[Abs[PhiAvAnalytical - PhiAvNumerical], 18]
];

Print["================================================"];