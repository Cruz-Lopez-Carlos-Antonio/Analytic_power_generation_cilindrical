(* ============================================================ *)
(* F0.wl                                                        *)
(* Operational implementation of F0(R)                          *)
(* ============================================================ *)

Module[{localDir = DirectoryName[$InputFileName]},
  Get[FileNameJoin[{localDir, "Parameters.wl"}]];
  Get[FileNameJoin[{localDir, "MR.wl"}]];
];

(* ============================================================ *)
(* PARAMETER                                                    *)
(* ============================================================ *)

PiD = alpha xi^2;

(* ============================================================ *)
(* NUMERICAL SETTINGS                                           *)
(* ============================================================ *)

wpF0 = 30;
agF0 = 14;
pgF0 = 14;

(* ============================================================ *)
(* AUXILIARY ODE                                                *)
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
    0 <= Rin <= 1, f0Sol[Rin],
    True, Indeterminate
];

(* ============================================================ *)
(* PUBLIC DERIVATIVE (OPERATIVA EXACTA)                         *)
(* ============================================================ *)

ClearAll[F0Prime];

F0Prime[Rin_?NumericQ] := Which[
    0 <= Rin <= 1, PiD Rin / (2 MR[Rin]),
    True, Indeterminate
];
