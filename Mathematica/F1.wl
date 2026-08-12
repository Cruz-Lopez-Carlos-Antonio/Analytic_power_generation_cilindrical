(* ============================================================ *)
(* F1.wl                                                        *)
(* Operational implementation of F1(R)                          *)
(* ============================================================ *)

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

epsQ = 10^-8;


(* ------------------------------------------------------------ *)
(* Auxiliary function                                           *)
(*                                                              *)
(* Q(R) = H(R)/(R M(R))                                         *)
(*                                                              *)
(* Near R = 0:                                                  *)
(* Q(R) = (Sinh[Psi(0)]/2) R + O(R^3)                          *)
(* ------------------------------------------------------------ *)

ClearAll[QF1];

QF1[Rin_?NumericQ] := Which[
    Rin == 0,
        0,

    0 < Rin < epsQ,
        (Sinh[centerPotential]/2) Rin,

    epsQ <= Rin <= 1,
        HR[Rin]/(Rin MR[Rin]),

    True,
        Indeterminate
];


(* ------------------------------------------------------------ *)
(* Auxiliary ODE                                                *)
(*                                                              *)
(* F1'(R) = Lambda delta^2 Q(R)                                 *)
(* F1(1)  = 0                                                   *)
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
    0 <= Rin <= 1,
        f1Sol[Rin],

    True,
        Indeterminate
];