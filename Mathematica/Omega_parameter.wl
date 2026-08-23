(* ============================================================ *)
(* Omega_parameter.wl                                           *)
(* Operational module for A, B and Omega                        *)
(* ============================================================ *)

omegaBaseDir = DirectoryName[$InputFileName];

Get[FileNameJoin[{omegaBaseDir, "F0.wl"}]];
Get[FileNameJoin[{omegaBaseDir, "F1.wl"}]];


(* ============================================================ *)
(* COMPATIBILITY ALIASES                                        *)
(* ============================================================ *)

If[Length[DownValues[F0]] == 0,
    F0[R_?NumericQ] := F0ODE[R]
];

If[Length[DownValues[F1]] == 0,
    F1[R_?NumericQ] := F1ODE[R]
];


(* ============================================================ *)
(* NUMERICAL SETTINGS                                           *)
(* Calibrated after validation of Omega_verification.wl         *)
(* ============================================================ *)

omegaWorkingPrecision = 30;
omegaAccuracyGoal = 12;
omegaPrecisionGoal = 12;
omegaMaxRecursion = 30;


(* ============================================================ *)
(* INTEGRANDS                                                   *)
(* ============================================================ *)

ClearAll[omegaIntegrandA, omegaIntegrandB];

omegaIntegrandA[R_?NumericQ] :=
    F0[R] Sinh[PsiPB[R]] R;

omegaIntegrandB[R_?NumericQ] :=
    F1[R] Sinh[PsiPB[R]] R;


(* ============================================================ *)
(* A AND B                                                      *)
(* ============================================================ *)

AOmega = NIntegrate[
    omegaIntegrandA[R],
    {R, 0, 1 - 1/delta, 1}, (* <-- División explícita en la capa límite física *)

    WorkingPrecision -> omegaWorkingPrecision,
    AccuracyGoal -> omegaAccuracyGoal,
    PrecisionGoal -> omegaPrecisionGoal,
    MaxRecursion -> omegaMaxRecursion,

    Method -> {
        "GlobalAdaptive",
        "SymbolicProcessing" -> 0
    }
];

BOmega = NIntegrate[
    omegaIntegrandB[R],
    {R, 0, 1 - 1/delta, 1}, (* <-- División explícita en la capa límite física *)

    WorkingPrecision -> omegaWorkingPrecision,
    AccuracyGoal -> omegaAccuracyGoal,
    PrecisionGoal -> omegaPrecisionGoal,
    MaxRecursion -> omegaMaxRecursion,

    Method -> {
        "GlobalAdaptive",
        "SymbolicProcessing" -> 0
    }
];


(* ============================================================ *)
(* OMEGA                                                        *)
(*                                                              *)
(* Omega = A/(1-B)                                              *)
(* ============================================================ *)

omegaDenominator = 1 - BOmega;

If[
    Abs[omegaDenominator] < 10^-8,
    Print[
        "WARNING in Omega_parameter.wl: 1 - B is numerically small: ",
        N[omegaDenominator, 18]
    ]
];

Omega = AOmega/omegaDenominator;
