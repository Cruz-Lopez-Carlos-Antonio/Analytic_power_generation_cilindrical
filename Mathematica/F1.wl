(* ============================================================ *)
(* F1.wl                                                        *)
(* Operational implementation of F1(R)                          *)
(* ============================================================ *)

Module[{localDir = DirectoryName[$InputFileName]},
  Get[FileNameJoin[{localDir, "HR.wl"}]];
  Get[FileNameJoin[{localDir, "MR.wl"}]];
  Get[FileNameJoin[{localDir, "Lambda_parameter.wl"}]];
];

(* ------------------------------------------------------------ *)
(* Numerical settings                                           *)
(* ------------------------------------------------------------ *)

wpF1 = 30;
agF1 = 14;
pgF1 = 14;
epsQ = 10^-8;

(* ------------------------------------------------------------ *)
(* Auxiliary function QF1(R)                                    *)
(* ------------------------------------------------------------ *)

ClearAll[QF1];

QF1[Rin_?NumericQ] := Which[
    Rin == 0, 0,
    0 < Rin < epsQ, (Sinh[centerPotential]/2) Rin,
    epsQ <= Rin <= 1, HR[Rin]/(Rin MR[Rin]),
    True, Indeterminate
];

(* ------------------------------------------------------------ *)
(* Auxiliary ODE                                                *)
(* ------------------------------------------------------------ *)

ClearAll[f1, f1Sol, F1];

f1Sol = NDSolveValue[
    {
        f1'[R] == Lambda delta^2 QF1[R],
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

(* ------------------------------------------------------------ *)
(* Public function                                              *)
(* ------------------------------------------------------------ *)

F1[Rin_?NumericQ] := Which[
    0 <= Rin <= 1, f1Sol[Rin],
    True, Indeterminate
];

(* ------------------------------------------------------------ *)
(* Public derivative (OPERATIVA EXACTA)                         *)
(* ------------------------------------------------------------ *)

ClearAll[F1Prime];

F1Prime[Rin_?NumericQ] := Which[
    0 <= Rin <= 1, Lambda delta^2 QF1[Rin],
    True, Indeterminate
];
