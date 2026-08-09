(* ========================================================= *)
(* Fcc.wl                                                     *)
(* ========================================================= *)
(*
   Computation and validation of

              F_cc = 2 Integral_0^1 R Cosh[Psi(R)] dR

   The previously validated Poisson-Boltzmann solution
   PsiPB[R] is loaded from:

              PoissonBoltzmann_.wl

   Two principal numerical methods are used:

       Method 1:
           Direct adaptive quadrature with NIntegrate.

       Method 2:
           Auxiliary differential equation

               Q'(R) = 2 R Cosh[Psi(R)],
               Q(0)  = 0,

           so that

               F_cc = Q(1).

   As an additional diagnostic, Method 1 is also repeated
   by splitting the integration interval.

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
      "PoissonBoltzmann_.wl"
    }
  ]
];


(* ========================================================= *)
(* 2. NUMERICAL SETTINGS                                     *)
(* ========================================================= *)

(*
   The quadrature accuracy goals are chosen conservatively
   relative to the precision of the previously computed
   Poisson-Boltzmann solution.

   WorkingPrecision controls internal arithmetic.
   It does not, by itself, imply that all 30 digits are
   physically or numerically reliable.
*)

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
(* 4. METHOD 1                                               *)
(*    DIRECT NUMERICAL QUADRATURE                            *)
(* ========================================================= *)

FccNIntegrate =
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


(* ========================================================= *)
(* 5. METHOD 2 - AUXILIARY ODE                                  *)
(* ========================================================= *)

ClearAll[Q];

(*
   Near the symmetry axis,

       Psi(R) = Psi(0) + O(R^2),

   hence

       Integral_0^epsPB 2 R Cosh[Psi(R)] dR
         = Cosh[Psi(0)] epsPB^2 + O(epsPB^4).

   This supplies the initial value for the cumulative
   integral without evaluating the PB InterpolatingFunction
   outside its domain.
*)

Qeps =
  Cosh[PsiPB[0]] epsPB^2;


QSolution =
  NDSolveValue[

    {
      Q'[R] ==
        2 R Cosh[PsiPB[R]],

      Q[epsPB] == Qeps
    },

    Q,

    {R, epsPB, 1},

    WorkingPrecision -> 40,

    AccuracyGoal -> 20,
    PrecisionGoal -> 20,

    MaxSteps -> Infinity
  ];


FccODE =
  QSolution[1];

(* ========================================================= *)
(* 6. ADDITIONAL VALIDATION                                  *)
(*    SPLIT NUMERICAL QUADRATURE                             *)
(* ========================================================= *)

(*
   Because Psi(R) varies most strongly near the wall,
   the interval is divided into two pieces.

   This is useful for checking whether the convergence
   warning from the direct quadrature is associated with
   adaptive refinement near R = 1.
*)

FccSplit1 =
  NIntegrate[

    fccIntegrand[R],

    {R, 0, 1/2},

    WorkingPrecision -> wpFcc,

    AccuracyGoal -> agFcc,
    PrecisionGoal -> pgFcc,

    Method -> {
      "GlobalAdaptive",
      "SymbolicProcessing" -> 0
    },

    MaxRecursion -> 50
  ];


FccSplit2 =
  NIntegrate[

    fccIntegrand[R],

    {R, 1/2, 1},

    WorkingPrecision -> wpFcc,

    AccuracyGoal -> agFcc,
    PrecisionGoal -> pgFcc,

    Method -> {
      "GlobalAdaptive",
      "SymbolicProcessing" -> 0
    },

    MaxRecursion -> 50
  ];


FccSplit =
  FccSplit1 + FccSplit2;


(* ========================================================= *)
(* 7. COMPARISON OF THE METHODS                              *)
(* ========================================================= *)

FccAbsDifference =
  Abs[
    FccNIntegrate - FccODE
  ];


FccRelDifference =
  FccAbsDifference/
    Max[
      Abs[FccNIntegrate],
      Abs[FccODE]
    ];


FccSplitAbsDifference =
  Abs[
    FccSplit - FccODE
  ];


FccSplitRelDifference =
  FccSplitAbsDifference/
    Max[
      Abs[FccSplit],
      Abs[FccODE]
    ];


(* ========================================================= *)
(* 8. DEFAULT VALUE FOR SUBSEQUENT MODULES                   *)
(* ========================================================= *)

(*
   For the moment, the direct quadrature is retained as
   the default Fcc value.

   After comparing Mathematica with the independent Python
   implementation, this choice can be reconsidered if
   necessary.
*)

Fcc =
  FccNIntegrate;


(* ========================================================= *)
(* 9. DIAGNOSTIC OUTPUT                                      *)
(* ========================================================= *)

Print[""];

Print[
  "=================================================="
];

Print[
  "F_cc validation"
];

Print[
  "=================================================="
];

Print[""];


Print[
  "Method 1 - Direct NIntegrate"
];

Print[
  "F_cc = ",
  N[FccNIntegrate, 18]
];

Print[""];


Print[
  "Method 2 - Auxiliary ODE"
];

Print[
  "F_cc = ",
  N[FccODE, 18]
];

Print[""];


Print[
  "Additional check - Split NIntegrate"
];

Print[
  "F_cc = ",
  N[FccSplit, 18]
];

Print[""];


Print[
  "--------------------------------------------------"
];

Print[
  "Direct quadrature vs ODE"
];

Print[
  "Absolute difference = ",
  ScientificForm[
    N[FccAbsDifference, 10]
  ]
];

Print[
  "Relative difference = ",
  ScientificForm[
    N[FccRelDifference, 10]
  ]
];

Print[""];


Print[
  "Split quadrature vs ODE"
];

Print[
  "Absolute difference = ",
  ScientificForm[
    N[FccSplitAbsDifference, 10]
  ]
];

Print[
  "Relative difference = ",
  ScientificForm[
    N[FccSplitRelDifference, 10]
  ]
];


Print[
  "--------------------------------------------------"
];

Print[""];

Print[
  "Default value used by subsequent modules:"
];

Print[
  "Fcc = ",
  N[Fcc, 18]
];

Print[""];

Print[
  "=================================================="
];