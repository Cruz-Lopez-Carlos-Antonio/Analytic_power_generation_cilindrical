(* ============================================================ *)
(* Omega_verification.wl                                        *)
(* Validation of A, B and Omega                                 *)
(* ============================================================ *)

baseDir = DirectoryName[$InputFileName];

Get[FileNameJoin[{baseDir, "F0.wl"}]];
Get[FileNameJoin[{baseDir, "F1.wl"}]];

(* ============================================================ *)
(* COMPATIBILITY ALIASES                                        *)
(* ============================================================ *)
(* Puente para asegurar que F0 y F1 estén siempre enlazados *)
If[Length[DownValues[F0]] == 0, F0[R_?NumericQ] := F0ODE[R]];
If[Length[DownValues[F1]] == 0, F1[R_?NumericQ] := F1ODE[R]];


(* ============================================================ *)
(* NUMERICAL SETTINGS                                           *)
(* ============================================================ *)

wpOmega = 30;
agOmega = 12; (* <-- Tolerancia calibrada para asimilar el ruido de interpolación *)
pgOmega = 12; (* <-- Tolerancia calibrada para asimilar el ruido de interpolación *)


(* ============================================================ *)
(* INTEGRANDS FOR A AND B                                       *)
(* ============================================================ *)

ClearAll[integrandA, integrandB];

integrandA[R_?NumericQ] :=
    F0[R] Sinh[PsiPB[R]] R;

integrandB[R_?NumericQ] :=
    F1[R] Sinh[PsiPB[R]] R;


(* ============================================================ *)
(* COMPUTE A AND B                                              *)
(*                                                              *)
(* A = Integral_0^1 F0(R) Sinh[Psi(R)] R dR                    *)
(* B = Integral_0^1 F1(R) Sinh[Psi(R)] R dR                    *)
(* ============================================================ *)

A = NIntegrate[
    integrandA[R],
    {R, 0, 1},

    WorkingPrecision -> wpOmega,
    AccuracyGoal -> agOmega,
    PrecisionGoal -> pgOmega,
    MaxRecursion -> 30, (* <-- El salvavidas numérico para evitar el límite de 9 *)

    Method -> {
        "GlobalAdaptive",
        "SymbolicProcessing" -> 0
    }
];

B = NIntegrate[
    integrandB[R],
    {R, 0, 1},

    WorkingPrecision -> wpOmega,
    AccuracyGoal -> agOmega,
    PrecisionGoal -> pgOmega,
    MaxRecursion -> 30, (* <-- El salvavidas numérico para evitar el límite de 9 *)

    Method -> {
        "GlobalAdaptive",
        "SymbolicProcessing" -> 0
    }
];


(* ============================================================ *)
(* COMPUTE OMEGA                                                *)
(*                                                              *)
(* Omega = A/(1-B)                                              *)
(* ============================================================ *)

denominatorOmega = 1 - B;

If[
    Abs[denominatorOmega] < 10^-8,
    Print[
        "WARNING: 1 - B is numerically small: ",
        N[denominatorOmega, 18]
    ]
];

Omega = A/denominatorOmega;


(* ============================================================ *)
(* RECONSTRUCT F FOR THE INTRINSIC CHECK                        *)
(* ============================================================ *)

ClearAll[FVerification];

FVerification[R_?NumericQ] :=
    F0[R] + Omega F1[R];


(* ============================================================ *)
(* DIRECT CHECK OF THE DEFINING FUNCTIONAL                      *)
(*                                                              *)
(* OmegaCheck = Integral_0^1                                    *)
(*              F(R) Sinh[Psi(R)] R dR                          *)
(* ============================================================ *)

ClearAll[integrandOmegaCheck];

integrandOmegaCheck[R_?NumericQ] :=
    FVerification[R] Sinh[PsiPB[R]] R;

OmegaCheck = NIntegrate[
    integrandOmegaCheck[R],
    {R, 0, 1},

    WorkingPrecision -> wpOmega,
    AccuracyGoal -> agOmega,
    PrecisionGoal -> pgOmega,
    MaxRecursion -> 30, (* <-- El salvavidas numérico para evitar el límite de 9 *)

    Method -> {
        "GlobalAdaptive",
        "SymbolicProcessing" -> 0
    }
];


(* ============================================================ *)
(* ERRORS                                                       *)
(* ============================================================ *)

OmegaAbsoluteError =
    Abs[OmegaCheck - Omega];

OmegaAPE =
    If[
        Abs[Omega] > 0,
        100 OmegaAbsoluteError/Abs[Omega],
        Indeterminate
    ];


(* ============================================================ *)
(* REPORT                                                       *)
(* ============================================================ *)

Print["=============================================="];
Print["A, B AND OMEGA VERIFICATION"];
Print["=============================================="];

Print["A              = ", N[A, 18]];
Print["B              = ", N[B, 18]];
Print["1 - B          = ", N[denominatorOmega, 18]];

Print[""];

Print["Omega          = ", N[Omega, 18]];
Print["OmegaCheck     = ", N[OmegaCheck, 18]];

Print[""];

Print["Absolute error = ", ScientificForm[OmegaAbsoluteError, 8]];
Print["APE (%)        = ", ScientificForm[OmegaAPE, 8]];

Print["=============================================="];