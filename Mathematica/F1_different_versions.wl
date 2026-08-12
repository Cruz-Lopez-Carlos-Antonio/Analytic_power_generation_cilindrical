(* ============================================================ *)
(* F1_validation.wl                                             *)
(* Validation of F1(R): direct quadrature vs auxiliary ODE      *)
(* ============================================================ *)

ClearAll["Global`*"];

baseDir = DirectoryName[$InputFileName];

Get[FileNameJoin[{baseDir, "HR.wl"}]];
Get[FileNameJoin[{baseDir, "MR.wl"}]];
Get[FileNameJoin[{baseDir, "Lambda_parameter.wl"}]];


(* ------------------------------------------------------------ *)
(* Numerical settings                                           *)
(* ------------------------------------------------------------ *)

wpF1 = 30;
agF1 = 15;
pgF1 = 15;

(* ------------------------------------------------------------ *)
(* Auxiliary function                                           *)
(*                                                              *)
(* Q(R) = H(R)/(R M(R))                                         *)
(*                                                              *)
(* The apparent singularity at R = 0 is removable, and          *)
(* Q(0) = 0.                                                    *)
(* ------------------------------------------------------------ *)

Q[Rin_?NumericQ] := Which[
    Rin == 0,
        0,
        
    0 < Rin < 10^-8,
        (Rin * Sinh[centerPotential]) / 2,

    10^-8 <= Rin <= 1,
        HR[Rin]/(Rin MR[Rin]),

    True,
        Indeterminate
];

(* ============================================================ *)
(* METHOD 1: DIRECT QUADRATURE                                   *)
(* ============================================================ *)

ClearAll[F1Direct];

F1Direct[Rin_?NumericQ] := Which[
    Rin == 1,
        0,

    0 <= Rin < 1,
        -Lambda delta^2 NIntegrate[
            Q[s],
            {s, Rin, 1},
            WorkingPrecision -> wpF1,
            AccuracyGoal -> agF1,
            PrecisionGoal -> pgF1,
            MaxRecursion -> 30, (* <-- Salvavidas numérico añadido *)
            Method -> {
                "GlobalAdaptive",
                "SymbolicProcessing" -> 0
            }
        ],

    True,
        Indeterminate
];

(* ============================================================ *)
(* METHOD 2: AUXILIARY ODE                                      *)
(*                                                              *)
(* F1'(R) = Lambda delta^2 Q(R)                                 *)
(* F1(1)  = 0                                                   *)
(* ============================================================ *)

ClearAll[f1, f1Sol, F1ODE];

f1Sol = NDSolveValue[
    {
        f1'[R] == Lambda delta^2 Q[R],
        f1[1] == 0
    },
    f1,
    {R, 0, 1},
    WorkingPrecision -> wpF1,
    AccuracyGoal -> agF1,
    PrecisionGoal -> pgF1,
    Method -> {
        "TimeIntegration" -> {
            "ExplicitRungeKutta",
            "DifferenceOrder" -> 8
        }
    }
];

F1ODE[Rin_?NumericQ] := Which[
    0 <= Rin <= 1,
        f1Sol[Rin],

    True,
        Indeterminate
];


(* ============================================================ *)
(* CONTROL POINTS                                               *)
(* ============================================================ *)

controlPoints = N[
    Table[k/19, {k, 0, 19}],
    wpF1
];


(* ------------------------------------------------------------ *)
(* Comparison table                                             *)
(* ------------------------------------------------------------ *)

comparisonData = Table[
    Module[
        {
            r,
            direct,
            ode,
            absError,
            ape
        },

        r = controlPoints[[k]];

        direct = F1Direct[r];
        ode = F1ODE[r];

        absError = Abs[direct - ode];

        ape = If[
            direct == 0,
            Missing["NotApplicable"],
            100 Abs[(ode - direct)/direct]
        ];

        {
            r,
            direct,
            ode,
            absError,
            ape
        }
    ],
    {k, Length[controlPoints]}
];


comparisonTable = Grid[
    Prepend[
        comparisonData,
        {
            "R",
            "F1Direct",
            "F1ODE",
            "Absolute error",
            "APE (%)"
        }
    ],
    Frame -> All,
    Alignment -> Right
];

Print[comparisonTable];


(* ============================================================ *)
(* ENDPOINT CHECKS                                              *)
(* ============================================================ *)

Print[""];
Print["--- Endpoint checks ---"];

Print[
    "F1Direct(1) = ",
    N[F1Direct[1], 18]
];

Print[
    "F1ODE(1)    = ",
    N[F1ODE[1], 18]
];

Print[
    "F1Direct(0) = ",
    N[F1Direct[0], 18]
];

Print[
    "F1ODE(0)    = ",
    N[F1ODE[0], 18]
];


(* ============================================================ *)
(* AXIS REGULARITY                                              *)
(*                                                              *)
(* Since Q(0)=0, F1'(0)=0.                                     *)
(* ============================================================ *)

Print[""];
Print["--- Axis regularity ---"];

Print[
    "Q(0) = ",
    Q[0]
];

Print[
    "Expected F1'(0) = Lambda delta^2 Q(0) = ",
    Lambda delta^2 Q[0]
];

Print[
    "ODE derivative at R=0 = ",
    N[f1Sol'[0], 18]
];


(* ============================================================ *)
(* DIFFERENTIAL RESIDUAL                                        *)
(*                                                              *)
(* residual(R) = F1_ODE'(R) - Lambda delta^2 Q(R)               *)
(* ============================================================ *)

ClearAll[F1Residual];

F1Residual[Rin_?NumericQ] :=
    f1Sol'[Rin] - Lambda delta^2 Q[Rin];


residualPoints = N[
    Table[k/20, {k, 1, 19}],
    wpF1
];

residualValues = Table[
    {
        r,
        F1Residual[r]
    },
    {r, residualPoints}
];

maxResidual = Max[
    Abs[residualValues[[All, 2]]]
];

Print[""];
Print["--- Differential residual ---"];

Print[
    "Maximum |F1'(R) - Lambda delta^2 Q(R)| = ",
    ScientificForm[N[maxResidual, 12]]
];


(* ============================================================ *)
(* GLOBAL COMPARISON                                            *)
(* ============================================================ *)

absoluteErrors =
    Cases[
        comparisonData[[All, 4]],
        _?NumericQ
    ];

apeValues =
    Cases[
        comparisonData[[All, 5]],
        _?NumericQ
    ];

Print[""];
Print["--- Global comparison ---"];

Print[
    "Maximum absolute error = ",
    ScientificForm[N[Max[absoluteErrors], 12]]
];

If[
    Length[apeValues] > 0,
    Print[
        "Maximum APE (%) = ",
        ScientificForm[N[Max[apeValues], 12]]
    ]
];


(* ============================================================ *)
(* PLOT                                                         *)
(* ============================================================ *)

Plot[
    {
        F1Direct[R],
        F1ODE[R]
    },
    {R, 0, 1},
    PlotLegends -> {
        "F1Direct",
        "F1ODE"
    },
    AxesLabel -> {
        "R",
        "F1(R)"
    },
    PlotRange -> All,
    PlotPoints -> 40,
    MaxRecursion -> 3,
    ImageSize -> Large
]