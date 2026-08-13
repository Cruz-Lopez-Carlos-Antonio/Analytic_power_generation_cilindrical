(* ============================================================ *)
(* F0_validation.wl                                             *)
(* Validation of F0(R): direct quadrature vs auxiliary ODE      *)
(* ============================================================ *)

ClearAll["Global`*"];

baseDir = DirectoryName[$InputFileName];

Get[FileNameJoin[{baseDir, "Parameters.wl"}]];
Get[FileNameJoin[{baseDir, "MR.wl"}]];


(* ============================================================ *)
(* PARAMETERS                                                   *)
(* ============================================================ *)

PiD = alpha xi^2;

wpF0 = 30;
agF0 = 12; (* <-- Exactitud calibrada para evitar slwcon *)
pgF0 = 12;


(* ============================================================ *)
(* METHOD 1: DIRECT QUADRATURE                                   *)
(* ============================================================ *)

ClearAll[F0Direct];

F0Direct[Rin_?NumericQ] := Which[

    Rin == 1,
        0,

    0 <= Rin < 1,
        -NIntegrate[
            PiD s/(2 MR[s]),
            {s, Rin, 1},

            WorkingPrecision -> wpF0,
            AccuracyGoal -> agF0,
            PrecisionGoal -> pgF0,

            Method -> {
                "GlobalAdaptive",
                "SymbolicProcessing" -> 0
            },

            MaxRecursion -> 30
        ],

    True,
        Indeterminate
];

(* ============================================================ *)
(* METHOD 2: AUXILIARY ODE                                      *)
(*                                                              *)
(* F0'(R) = PiD R/(2 M(R))                                      *)
(* F0(1)  = 0                                                   *)
(* ============================================================ *)

ClearAll[f0, f0Func, f0Prime, F0ODE];

f0Func = NDSolveValue[
    {
        f0'[R] == PiD R/(2 MR[R]),
        f0[1] == 0
    },
    f0, (* <-- Pedimos ÚNICAMENTE la función para evitar ceros simbólicos *)
    {R, 0, 1},

    WorkingPrecision -> wpF0,
    AccuracyGoal -> agF0,
    PrecisionGoal -> pgF0,

    Method -> {
        "TimeIntegration" -> {
            "ExplicitRungeKutta",
            "DifferenceOrder" -> 8
        }
    }
];

(* Extraemos la derivada de forma robusta y nativa *)
f0Prime = Derivative[1][f0Func];

F0ODE[Rin_?NumericQ] := Which[
    0 <= Rin <= 1,
        f0Func[Rin],
    True,
        Indeterminate
];

(* ============================================================ *)
(* CONTROL POINTS                                               *)
(* ============================================================ *)

controlPoints = N[
    Table[k/19, {k, 0, 19}],
    wpF0
];


(* ============================================================ *)
(* COMPARISON TABLE                                             *)
(* ============================================================ *)

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

        direct = F0Direct[r];
        ode = F0ODE[r];

        absError = Abs[
            direct - ode
        ];

        ape = If[
            direct == 0,
            Missing["NotApplicable"],
            100 Abs[
                (ode - direct)/direct
            ]
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
            "F0Direct",
            "F0ODE",
            "Absolute error",
            "APE (%)"
        }
    ],

    Frame -> All,
    Alignment -> Right
];


Print[
    comparisonTable
];


(* ============================================================ *)
(* ENDPOINT CHECKS                                              *)
(* ============================================================ *)

Print[""];
Print["--- Endpoint checks ---"];

Print[
    "PiD         = ",
    N[PiD, 18]
];

Print[
    "F0Direct(1) = ",
    N[F0Direct[1], 18]
];

Print[
    "F0ODE(1)    = ",
    N[F0ODE[1], 18]
];

Print[
    "F0Direct(0) = ",
    N[F0Direct[0], 18]
];

Print[
    "F0ODE(0)    = ",
    N[F0ODE[0], 18]
];


(* ============================================================ *)
(* AXIS REGULARITY                                              *)
(*                                                              *)
(* Since M(0)=1,                                                *)
(*                                                              *)
(* F0'(0) = PiD*0/(2 M(0)) = 0.                                *)
(* ============================================================ *)

Print[""];
Print["--- Axis regularity ---"];

Print[
    "M(0) = ",
    N[MR[0], 18]
];

Print[
    "Expected F0'(0) = ",
    N[PiD 0/(2 MR[0]), 18]
];

Print[
    "ODE derivative at R=0 = ",
    N[f0Prime[0], 18]
];


(* ============================================================ *)
(* DIFFERENTIAL RESIDUAL                                        *)
(* ============================================================ *)

ClearAll[F0Residual];

F0Residual[Rin_?NumericQ] :=
    f0Prime[Rin] - PiD Rin/(2 MR[Rin]);


residualPoints = N[
    Table[k/20, {k, 1, 19}],
    wpF0
];


residualValues = Table[
    {
        r,
        F0Residual[r]
    },
    {r, residualPoints}
];


maxResidual = Max[
    Abs[
        residualValues[[All, 2]]
    ]
];


Print[""];
Print["--- Differential residual ---"];

Print[
    "Maximum |F0'(R) - PiD R/(2 M(R))| = ",
    ScientificForm[
        N[maxResidual, 12]
    ]
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
    ScientificForm[
        N[
            Max[absoluteErrors],
            12
        ]
    ]
];


If[
    Length[apeValues] > 0,

    Print[
        "Maximum APE (%) = ",
        ScientificForm[
            N[
                Max[apeValues],
                12
            ]
        ]
    ]
];


(* ============================================================ *)
(* SIGN CHECK                                                   *)
(* ============================================================ *)

signValues = Table[
    F0ODE[r],
    {r, Most[controlPoints]}
];


Print[""];
Print["--- Sign check ---"];

Print[
    "F0(R) < 0 for 0 <= R < 1 : ",
    And @@ Thread[
        signValues < 0
    ]
];


(* ============================================================ *)
(* PLOT                                                         *)
(* ============================================================ *)

Plot[
    {
        F0Direct[R],
        F0ODE[R]
    },

    {R, 0, 1},

    PlotLegends -> {
        "F0Direct",
        "F0ODE"
    },

    AxesLabel -> {
        "R",
        "F0(R)"
    },

    PlotRange -> All,
    PlotPoints -> 40,
    MaxRecursion -> 3,
    ImageSize -> Large
]