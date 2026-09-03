(* ============================================================ *)
(* Code developed by Sánchez Lozano, G., Cruz-López C.-A., and  *)
(* F. Méndez, for the research:                                 *)
(* "Power generation in cylindrical microchannels with high     *)
(* surface zeta potential: a survey on viscoelectric effects    *)
(* over streaming potentials and efficiency"                    *)
(* Year: 2026                                                   *)
(* ------------------------------------------------------------ *)
(* F_cc_parametric.wl                                           *)
(* ============================================================ *)
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
agFccParam = 12; (* Harmonized tolerance for baseline noise *)
pgFccParam = 12; (* Harmonized tolerance for baseline noise *)


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
  deltaLoc = pb["Delta"]; (* Dynamic extraction of the physical parameter *)


  integrand[R_?NumericQ] :=
   2 R Cosh[psi[R]];


  fccValue =
   NIntegrate[
    integrand[R],

    {R, 0, 1 - 1/deltaLoc, 1}, (* Explicit mapping of the boundary layer *)

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
