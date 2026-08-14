(* ============================================================ *)
(* G.wl                                                        *)
(* Dimensionless modified-pressure gradient                    *)
(* ============================================================ *)

Module[{localDir = DirectoryName[$InputFileName]},
  Get[FileNameJoin[{localDir, "F.wl"}]];
];

(* ============================================================ *)
(* NUMERICAL SETTINGS                                           *)
(* ============================================================ *)

wpG = 30;
agG = 12;
pgG = 12;

(* ============================================================ *)
(* FLOW INTEGRAL                                                *)
(* ============================================================ *)

ClearAll[IF];

IF = NIntegrate[
    R F[R],
    {R, 0, 1},
    WorkingPrecision -> wpG,
    AccuracyGoal -> agG,
    PrecisionGoal -> pgG,
    MaxRecursion -> 30,
    Method -> "GlobalAdaptive"
];

(* ============================================================ *)
(* DIMENSIONLESS PRESSURE GRADIENT                              *)
(* ============================================================ *)

ClearAll[G];

G = 1/(2 IF);