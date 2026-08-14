(* ============================================================ *)
(* G_verification.wl                                            *)
(* Verification of the dimensionless pressure gradient G        *)
(* ============================================================ *)

Module[{localDir = DirectoryName[$InputFileName]},
  Get[FileNameJoin[{localDir, "F.wl"}]];
];

(* ============================================================ *)
(* NUMERICAL SETTINGS                                           *)
(* ============================================================ *)

wpG = 30;
agG = 12; (* Calibración armónica heredada de F(R) *)
pgG = 12; (* Calibración armónica heredada de F(R) *)
maxRecursionG = 30; (* Salvavidas para evitar el límite de 9 *)

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
    MaxRecursion -> maxRecursionG,
    Method -> "GlobalAdaptive"
];

(* ============================================================ *)
(* NON-DEGENERACY CHECK                                         *)
(* ============================================================ *)

If[
    !NumericQ[IF] ||
    MemberQ[{Indeterminate, ComplexInfinity, Infinity, -Infinity}, IF],
    Print["ERROR: the integral IF is not a finite numerical value."];
    Abort[];
];

If[
    Abs[IF] < 10^-12,
    Print["ERROR: the denominator is numerically degenerate."];
    Abort[];
];

(* ============================================================ *)
(* PRESSURE GRADIENT                                            *)
(* ============================================================ *)

ClearAll[G];

G = 1/(2 IF);

(* ============================================================ *)
(* FLOW-RATE NORMALIZATION CHECK                                *)
(* ============================================================ *)

ClearAll[flowCheck, flowResidual];

flowCheck = NIntegrate[
    R G F[R],
    {R, 0, 1},
    WorkingPrecision -> wpG,
    AccuracyGoal -> agG,
    PrecisionGoal -> pgG,
    MaxRecursion -> maxRecursionG,
    Method -> "GlobalAdaptive"
];

flowResidual = flowCheck - 1/2;

(* ============================================================ *)
(* OUTPUT                                                       *)
(* ============================================================ *)

Print["=============================================="];
Print["Verification of pressure gradient G"];
Print["=============================================="];

Print["Integral IF = Integrate_0^1 R F(R) dR"];
Print["IF = ", N[IF, 18]];

Print[""];
Print["G = 1/(2 IF)"];
Print["G = ", N[G, 18]];

Print[""];
Print["Flow-rate normalization"];
Print["Integral_0^1 R G F(R) dR = ", N[flowCheck, 18]];
Print["Target                         = ", N[1/2, 18]];
Print["Absolute residual              = ", ScientificForm[Abs[flowResidual], 8]];

Print[""];
Print["Denominator magnitude |IF| = ", ScientificForm[Abs[IF], 8]];
Print["Sign[IF] = ", Sign[IF]];
Print["Sign[G]  = ", Sign[G]];

Print["=============================================="];