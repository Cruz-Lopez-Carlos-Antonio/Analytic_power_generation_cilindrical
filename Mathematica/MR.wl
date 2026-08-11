(* ========================================================= *)
(* M(R) — VISCOELECTRIC FACTOR                               *)
(* ========================================================= *)
(*
   Defines the dimensionless function

       M(R) = Exp[omega (Psi'(R))^2],

   where omega is obtained from Parameters.wl and Psi'(R)
   from the validated Poisson-Boltzmann solution.

   Since Psi'(0) = 0 by symmetry,

       M(0) = 1.

   This module is intended to remain silent when loaded.
*)
(* ========================================================= *)


(* ========================================================= *)
(* 1. LOAD POISSON-BOLTZMANN SOLUTION                        *)
(* ========================================================= *)

moduleDirectory = DirectoryName[$InputFileName];

Get[
  FileNameJoin[
    {moduleDirectory, "PoissonBoltzmann.wl"}
  ]
];


(* ========================================================= *)
(* 2. DEFINE M(R)                                            *)
(* ========================================================= *)

ClearAll[MR];

MR[Rin_?NumericQ] :=
  Which[
    Rin == 0,
    1,

    0 < Rin <= 1,
    Exp[
      omega*(PsiPrimePB[Rin])^2
    ],

    True,
    Indeterminate
  ];
