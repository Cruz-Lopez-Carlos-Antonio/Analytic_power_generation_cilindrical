(* ========================================================= *)
(* H(R) — AUXILIARY ODE                                      *)
(* ========================================================= *)

moduleDirectory =
  DirectoryName[$InputFileName];

Get[
  FileNameJoin[
    {
      moduleDirectory,
      "PoissonBoltzman.wl"
    }
  ]
];


(* ========================================================= *)
(* AUXILIARY ODE                                             *)
(* ========================================================= *)
(*
   H(R) = Integral_0^R tau Sinh[Psi(tau)] d tau

   Therefore:

       H'(R) = R Sinh[Psi(R)],
       H(0)  = 0.
*)
(* ========================================================= *)

ClearAll[h];

hSolution =
  NDSolveValue[
    {
      h'[R] ==
        R Sinh[PsiPB[R]],

      h[0] == 0
    },

    h,

    {
      R,
      0,
      1
    },

    WorkingPrecision -> 30,
    AccuracyGoal -> 18,
    PrecisionGoal -> 18,

    MaxSteps -> Infinity
  ];


(* ========================================================= *)
(* PUBLIC FUNCTION H(R)                                      *)
(* ========================================================= *)

ClearAll[HR];

HR[Rin_?NumericQ] :=
  Which[

    Rin == 0,
    0,

    0 < Rin <= 1,
    hSolution[Rin],

    True,
    Indeterminate
  ];