(* ========================================================= *)
(* F_cc_parametric.wl                                        *)
(* ========================================================= *)
(*
   Parametric computation of

       F_cc = 2 Integral_0^1 R Cosh[Psi(R)] dR

   using a previously constructed Poisson-Boltzmann
   solution object:

       pb = SolvePB[deltaValue, psiSValue];

   The Poisson-Boltzmann problem is NOT solved again here.
*)
(* ========================================================= *)


(* ========================================================= *)
(* 1. NUMERICAL SETTINGS                                     *)
(* ========================================================= *)

wpFccParam = 30;
agFccParam = 12; (* Tolerancia armonizada para ruido basal *)
pgFccParam = 12; (* Tolerancia armonizada para ruido basal *)


(* ========================================================= *)
(* 2. PARAMETRIC F_cc                                        *)
(* ========================================================= *)

ClearAll[ComputeFcc];

ComputeFcc[pb_Association] :=
 Module[
  {
   psi,
   deltaLoc,
   integrand,
   fccValue
   },

  psi = pb["Psi"];
  deltaLoc = pb["Delta"]; (* Extracción dinámica del parámetro físico *)


  integrand[R_?NumericQ] :=
   2 R Cosh[psi[R]];


  fccValue =
   NIntegrate[
    integrand[R],

    {R, 0, 1 - 1/deltaLoc, 1}, (* Mapeo explícito de la capa límite *)

    WorkingPrecision -> wpFccParam,

    AccuracyGoal -> agFccParam,
    PrecisionGoal -> pgFccParam,

    Method -> {
      "GlobalAdaptive",
      "SymbolicProcessing" -> 0
      },

    MaxRecursion -> 50
    ];


  fccValue
  ];
