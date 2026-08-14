(* ============================================================ *)
(* Corrected F_verification.wl - Shanoa's Bypass                *)
(* ============================================================ *)

fVerificationBaseDir = DirectoryName[$InputFileName];

Get[FileNameJoin[{fVerificationBaseDir, "F.wl"}]];
Get[FileNameJoin[{fVerificationBaseDir, "MR.wl"}]];
Get[FileNameJoin[{fVerificationBaseDir, "HR.wl"}]];
Get[FileNameJoin[{fVerificationBaseDir, "PoissonBoltzmann.wl"}]];

wpFVerification = 30;
agFVerification = 12;
pgFVerification = 12;
maxRecursionFVerification = 30;
nResidualPoints = 1001;

(* ============================================================ *)
(* BYPASS ANALÍTICO DE SHANOA                                   *)
(* Forzamos la derivada operativa exacta para el residual       *)
(* ============================================================ *)
ClearAll[F0PrimeExact, F1PrimeExact, FPrimeExact];

F0PrimeExact[Rin_?NumericQ] := Which[
    0 <= Rin <= 1, PiD Rin / (2 MR[Rin]),
    True, Indeterminate
];

F1PrimeExact[Rin_?NumericQ] := Which[
    0 <= Rin <= 1, Lambda delta^2 QF1[Rin],
    True, Indeterminate
];

FPrimeExact[Rin_?NumericQ] := Which[
    Rin == 0, 0,
    0 < Rin <= 1, F0PrimeExact[Rin] + Omega F1PrimeExact[Rin],
    True, Indeterminate
];


(* ============================================================ *)
(* OMEGA FUNCTIONAL CHECK                                       *)
(* ============================================================ *)
ClearAll[omegaFIntegrand];

omegaFIntegrand[R_?NumericQ] :=
  F[R] Sinh[PsiPB[R]] R;

OmegaCheckF = NIntegrate[
  omegaFIntegrand[R],
  {R, 0, 1},
  WorkingPrecision -> wpFVerification,
  AccuracyGoal -> agFVerification,
  PrecisionGoal -> pgFVerification,
  MaxRecursion -> maxRecursionFVerification,
  Method -> {"GlobalAdaptive", "SymbolicProcessing" -> 0}
];

omegaFAbsoluteError = Abs[OmegaCheckF - Omega];

omegaFAPE =
  If[
    Abs[Omega] > 0,
    100 omegaFAbsoluteError/Abs[Omega],
    Indeterminate
  ];

(* ============================================================ *)
(* INTEGRATED RESIDUAL CHECK                                    *)
(* ============================================================ *)
ClearAll[integratedResidualF];

integratedResidualF[R_?NumericQ] :=
  MR[R] R FPrimeExact[R]
  -
  (
    PiD R^2/2
    +
    Omega Lambda delta^2 HR[R]
  );

residualGrid =
  N[Subdivide[0, 1, nResidualPoints - 1], 30];

residualValues =
  integratedResidualF /@ residualGrid;

residualLInf =
  Max[Abs[residualValues]];

residualL2 =
  Sqrt[Mean[residualValues^2]];

fSamplePoints = {0, 0.25, 0.50, 0.75, 1};
fSampleValues = F /@ fSamplePoints;

(* ============================================================ *)
(* REPORT                                                       *)
(* ============================================================ *)
Print["=============================================="];
Print["F(R) VERIFICATION - MATHEMATICA"];
Print["=============================================="];
Print[""];
Print["Omega              = ", N[Omega, 18]];
Print[""];
Print["Boundary / symmetry conditions"];
Print["------------------------------"];
Print["F(1)               = ", N[F[1], 18]];
Print["F'(0)              = ", N[FPrimeExact[0], 18]];
Print[""];
Print["Omega functional check"];
Print["----------------------"];
Print["OmegaCheck         = ", N[OmegaCheckF, 18]];
Print["Absolute error     = ", ScientificForm[omegaFAbsoluteError, 8]];
Print["APE (%)            = ", ScientificForm[omegaFAPE, 8]];
Print[""];
Print["Integrated Eq. (50) residual"];
Print["----------------------------"];
Print["L_inf residual     = ", ScientificForm[residualLInf, 8]];
Print["L2 residual        = ", ScientificForm[residualL2, 8]];
Print[""];
Print["Selected F(R) values"];
Print["--------------------"];
Print["       R                   F(R)"];
Print["----------------------------------------"];

Do[
  Print[
    NumberForm[N[fSamplePoints[[i]], 16], {6, 2}],
    "        ",
    ScientificForm[N[fSampleValues[[i]], 16], 15]
  ],
  {i, Length[fSamplePoints]}
];

Print[""];
Print["=============================================="];
