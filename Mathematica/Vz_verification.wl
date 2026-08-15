(* ============================================================ *)
(* Vz_verification.wl                                          *)
(* Verification of the dimensionless axial velocity V_Z(R)     *)
(* ============================================================ *)

Module[{localDir = DirectoryName[$InputFileName]},
  Get[FileNameJoin[{localDir, "G.wl"}]];
];

(* ============================================================ *)
(* NUMERICAL SETTINGS                                           *)
(* ============================================================ *)

wpVz = 30;
agVz = 12;
pgVz = 12;

(* ============================================================ *)
(* AXIAL VELOCITY                                               *)
(* ============================================================ *)

ClearAll[Vz, VzPrime];

Vz[R_?NumericQ] := G F[R];

(* Operational derivative:
 *
 *     Vz'(R) = G F'(R)
 *
 * We use the already validated operational derivative FPrime[R]
 * instead of numerically differentiating F(R).
 *)

VzPrime[R_?NumericQ] := G FPrime[R];

(* ============================================================ *)
(* BOUNDARY AND REGULARITY CHECKS                               *)
(* ============================================================ *)

wallValue = Vz[1];

axisDerivative = VzPrime[0];

wallResidual = Abs[wallValue];

axisResidual = Abs[axisDerivative];

(* ============================================================ *)
(* FLOW-RATE NORMALIZATION                                      *)
(* ============================================================ *)

flowCheck = NIntegrate[
    R Vz[R],
    {R, 0, 1},
    WorkingPrecision -> wpVz,
    AccuracyGoal -> agVz,
    PrecisionGoal -> pgVz,
    MaxRecursion -> 30,
    Method -> "GlobalAdaptive"
];

flowTarget = 1/2;

flowResidual = Abs[
    flowCheck - flowTarget
];

(* ============================================================ *)
(* CONSISTENCY CHECK USING G AND IF                             *)
(* ============================================================ *)

flowFromDefinition = G IF;

definitionResidual = Abs[
    flowFromDefinition - 1/2
];

(* ============================================================ *)
(* REFERENCE PROFILE POINTS                                     *)
(* ============================================================ *)

samplePoints = {
    0,
    1/4,
    1/2,
    3/4,
    1
};

sampleTable = Table[
    {
        R,
        Vz[N[R, wpVz]]
    },
    {R, samplePoints}
];

(* ============================================================ *)
(* OUTPUT                                                       *)
(* ============================================================ *)

Print["=============================================="];
Print["Verification of axial velocity V_Z(R)"];
Print["=============================================="];

Print[""];
Print["Pressure gradient"];
Print["G = ", N[G, 18]];

Print[""];
Print["Boundary and regularity conditions"];
Print["V_Z(1)     = ", N[wallValue, 18]];
Print["|V_Z(1)|   = ", ScientificForm[wallResidual, 8]];

Print[""];
Print["V_Z'(0)    = ", N[axisDerivative, 18]];
Print["|V_Z'(0)|  = ", ScientificForm[axisResidual, 8]];

Print[""];
Print["Flow-rate normalization"];
Print[
    "Integral_0^1 R V_Z(R) dR = ",
    N[flowCheck, 18]
];
Print[
    "Target                    = ",
    N[flowTarget, 18]
];
Print[
    "Absolute residual         = ",
    ScientificForm[flowResidual, 8]
];

Print[""];
Print["Consistency from G IF"];
Print[
    "G IF                      = ",
    N[flowFromDefinition, 18]
];
Print[
    "Absolute residual         = ",
    ScientificForm[definitionResidual, 8]
];

Print[""];
Print["Reference profile"];
Print[
    Grid[
        Prepend[
            sampleTable,
            {"R", "V_Z(R)"}
        ],
        Frame -> All
    ]
];