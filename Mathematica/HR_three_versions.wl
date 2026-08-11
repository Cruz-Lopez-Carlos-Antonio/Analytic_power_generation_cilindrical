(* ========================================================= *)
(* H(R) — DIFFERENT NUMERICAL VERSIONS                       *)
(* ========================================================= *)
(*
   We define

       H(R) = Integral_0^R tau Sinh[Psi(tau)] d tau

   and compare three numerical approaches:

   1. Direct NIntegrate
   2. Auxiliary ODE
   3. Cumulative trapezoidal integration + interpolation

   This file is intended for numerical validation and
   comparison of the three approaches.
*)
(* ========================================================= *)


(* ========================================================= *)
(* 1. LOAD POISSON-BOLTZMANN SOLUTION                        *)
(* ========================================================= *)

moduleDirectory =
  DirectoryName[$InputFileName];

Get[
  FileNameJoin[
    {
      moduleDirectory,
      "PoissonBoltzmann.wl"
    }
  ]
];


(* ========================================================= *)
(* 2. SAFE EVALUATION OF Psi(R)                              *)
(* ========================================================= *)
(*
   Some numerical algorithms may internally evaluate the
   interpolating function at points slightly outside [0,1].

   PsiSafe prevents accidental extrapolation by projecting
   such numerical evaluation points back onto the physical
   interval.
*)
(* ========================================================= *)

ClearAll[PsiSafe];

PsiSafe[x_?NumericQ] :=
  PsiPB[
    Min[
      1,
      Max[
        0,
        x
      ]
    ]
  ];


(* ========================================================= *)
(* 3. METHOD 1 — DIRECT NINTEGRATE                           *)
(* ========================================================= *)

ClearAll[HDirect];

ClearAll[HDirect];

HDirect[Rin_?NumericQ] :=
  Which[

    Rin == 0,
    0,

    0 < Rin <= 1,
    NIntegrate[
      tau Sinh[PsiSafe[tau]],
      {
        tau,
        0,
        Rin
      },

      WorkingPrecision -> 30,
      AccuracyGoal -> 18,
      PrecisionGoal -> 18,
      MaxRecursion -> 30, (* <-- Bifurcations' limit *)
      Method -> "GlobalAdaptive"
    ],

    True,
    Indeterminate
  ];


(* ========================================================= *)
(* 4. METHOD 2 — AUXILIARY ODE                               *)
(* ========================================================= *)
(*
   Since

       H(R) = Integral_0^R tau Sinh[Psi(tau)] d tau,

   we have

       H'(R) = R Sinh[Psi(R)],
       H(0)  = 0.

   Therefore H can be obtained from a first-order IVP.
*)
(* ========================================================= *)

ClearAll[h];

hSolution =
  NDSolveValue[
    {
      h'[R] ==
        R Sinh[PsiSafe[R]],

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


ClearAll[HODE];

HODE[Rin_?NumericQ] :=
  Which[

    Rin == 0,
    0,

    0 < Rin <= 1,
    hSolution[Rin],

    True,
    Indeterminate
  ];


(* ========================================================= *)
(* 5. METHOD 3 — CUMULATIVE TRAPEZOID + INTERPOLATION        *)
(* ========================================================= *)

nGrid = 5000;


rGrid =
  N[
    Subdivide[
      0,
      1,
      nGrid
    ],
    30
  ];


gGrid =
  Table[
    r Sinh[PsiSafe[r]],
    {
      r,
      rGrid
    }
  ];


hGrid =
  Prepend[
    Accumulate[
      Table[

        (
          rGrid[[i + 1]]
          -
          rGrid[[i]]
        )
        *
        (
          gGrid[[i + 1]]
          +
          gGrid[[i]]
        )
        /2,

        {
          i,
          1,
          Length[rGrid] - 1
        }
      ]
    ],

    0
  ];


hInterpolation =
  Interpolation[
    Transpose[
      {
        rGrid,
        hGrid
      }
    ],

    InterpolationOrder -> 3
  ];


ClearAll[HTrap];

HTrap[Rin_?NumericQ] :=
  Which[

    Rin == 0,
    0,

    0 < Rin <= 1,
    hInterpolation[Rin],

    True,
    Indeterminate
  ];


(* ========================================================= *)
(* END OF DEFINITIONS                                        *)
(* ========================================================= *)
