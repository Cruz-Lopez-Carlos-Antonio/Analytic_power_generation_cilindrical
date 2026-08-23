(* ============================================================ *)
(* Omega_parameter_parametric.wl                                *)
(* Parametric operational module for A, B and Omega             *)
(* ============================================================ *)

ClearAll[SolveOmegaParametric];

SolveOmegaParametric[
   pb_Association,
   f0Obj_Association,
   f1Obj_Association,
   deltaInput_?NumericQ
   ] :=
 Module[
  {
   omegaIntegrandAParametric,
   omegaIntegrandBParametric,
   AParametric,
   BParametric,
   denominatorParametric,
   OmegaParametric,
   rSplit,
   integrationRange
   },

  omegaWorkingPrecisionParametric = 30;
  omegaAccuracyGoalParametric = 12;
  omegaPrecisionGoalParametric = 12;
  omegaMaxRecursionParametric = 30;

  If[
   deltaInput > 1,
   rSplit = 1 - 1/deltaInput;
   integrationRange = {R, 0, rSplit, 1},
   rSplit = None;
   integrationRange = {R, 0, 1}
   ];

  ClearAll[
   omegaIntegrandAParametric,
   omegaIntegrandBParametric
   ];

  omegaIntegrandAParametric[R_?NumericQ] :=
   f0Obj["F0"][R]*
    Sinh[pb["Psi"][R]]*
    R;

  omegaIntegrandBParametric[R_?NumericQ] :=
   f1Obj["F1"][R]*
    Sinh[pb["Psi"][R]]*
    R;

  AParametric =
   NIntegrate[
    omegaIntegrandAParametric[R],
    Evaluate[integrationRange],
    WorkingPrecision -> omegaWorkingPrecisionParametric,
    AccuracyGoal -> omegaAccuracyGoalParametric,
    PrecisionGoal -> omegaPrecisionGoalParametric,
    MaxRecursion -> omegaMaxRecursionParametric,
    Method -> {
      "GlobalAdaptive",
      "SymbolicProcessing" -> 0
      }
    ];

  BParametric =
   NIntegrate[
    omegaIntegrandBParametric[R],
    Evaluate[integrationRange],
    WorkingPrecision -> omegaWorkingPrecisionParametric,
    AccuracyGoal -> omegaAccuracyGoalParametric,
    PrecisionGoal -> omegaPrecisionGoalParametric,
    MaxRecursion -> omegaMaxRecursionParametric,
    Method -> {
      "GlobalAdaptive",
      "SymbolicProcessing" -> 0
      }
    ];

  denominatorParametric = 1 - BParametric;

  If[
   Abs[denominatorParametric] < 10^-8,
   Print[
    "WARNING in Omega_parameter_parametric.wl: ",
    "1 - B is numerically small: ",
    N[denominatorParametric, 18]
    ]
   ];

  OmegaParametric =
   AParametric/denominatorParametric;

  <|
   "A" -> AParametric,
   "B" -> BParametric,
   "Denominator" -> denominatorParametric,
   "Omega" -> OmegaParametric,
   "Delta" -> deltaInput,
   "SplitPoint" -> rSplit
   |>
  ];
