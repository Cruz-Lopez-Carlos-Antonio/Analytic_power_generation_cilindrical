(* ============================================================ *)
(* F0.wl                                                        *)
(* Operational implementation of F0(R)                          *)
(* ============================================================ *)

baseDir = DirectoryName[$InputFileName];

Get[FileNameJoin[{baseDir, "Parameters.wl"}]];
Get[FileNameJoin[{baseDir, "MR.wl"}]];


(* ============================================================ *)
(* PARAMETER                                                    *)
(* ============================================================ *)

PiD = alpha xi^2;


(* ============================================================ *)
(* NUMERICAL SETTINGS                                           *)
(* ============================================================ *)

wpF0 = 30;
agF0 = 15;
pgF0 = 15;


(* ============================================================ *)
(* AUXILIARY ODE                                                *)
(*                                                              *)
(* F0'(R) = PiD R/(2 M(R))                                      *)
(* F0(1)  = 0                                                   *)
(* ============================================================ *)

ClearAll[f0, f0Sol, F0];

f0Sol = NDSolveValue[
    {
        f0'[R] == PiD R/(2 MR[R]),
        f0[1] == 0
    },
    f0,
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


(* ============================================================ *)
(* PUBLIC FUNCTION                                              *)
(* ============================================================ *)

F0[Rin_?NumericQ] := Which[
    0 <= Rin <= 1,
        f0Sol[Rin],

    True,
        Indeterminate
];


(* ============================================================ *)
(* OPTIONAL PUBLIC DERIVATIVE                                   *)
(* ============================================================ *)

f0Prime = Derivative[1][f0Sol];

F0Prime[Rin_?NumericQ] := Which[
    0 <= Rin <= 1,
        f0Prime[Rin],

    True,
        Indeterminate
];