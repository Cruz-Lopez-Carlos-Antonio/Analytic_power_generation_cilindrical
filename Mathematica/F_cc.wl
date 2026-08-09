(* ========================================================= *)
(* Fcc.wl                                                     *)
(* ========================================================= *)
(*
   Computation of

       F_cc = 2 Integral_0^1 R Cosh[Psi(R)] dR

   using the previously validated Poisson-Boltzmann
   solution PsiPB[R].
*)
(* ========================================================= *)


(* ========================================================= *)
(* 1. LOAD POISSON-BOLTZMANN MODULE                          *)
(* ========================================================= *)

moduleDirectory = DirectoryName[$InputFileName];

Get[
  FileNameJoin[
    {
      moduleDirectory,
      "PoissonBoltzmann.wl"
    }
  ]
];


(* ========================================================= *)
(* 2. NUMERICAL SETTINGS                                     *)
(* ========================================================= *)

wpFcc = 30;
agFcc = 14;
pgFcc = 14;


(* ========================================================= *)
(* 3. INTEGRAND                                              *)
(* ========================================================= *)

ClearAll[fccIntegrand];

fccIntegrand[R_?NumericQ] :=
  2 R Cosh[PsiPB[R]];


(* ========================================================= *)
(* 4. COMPUTE F_cc                                           *)
(* ========================================================= *)

Fcc =
  NIntegrate[

    fccIntegrand[R],

    {R, 0, 1},

    WorkingPrecision -> wpFcc,

    AccuracyGoal -> agFcc,
    PrecisionGoal -> pgFcc,

    Method -> {
      "GlobalAdaptive",
      "SymbolicProcessing" -> 0
    },

    MaxRecursion -> 50
  ];

