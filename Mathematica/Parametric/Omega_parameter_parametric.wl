(* ============================================================ *)
(* Code developed by Sánchez Lozano, G., Cruz-López C.-A., and  *)
(* F. Méndez, for the research:                                 *)
(* "Power generation in cylindrical microchannels with high     *)
(* surface zeta potential: a survey on viscoelectric effects    *)
(* over streaming potentials and efficiency"                    *)
(* Year: 2026                                                   *)
(* ------------------------------------------------------------ *)
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
   rangeInner,
   rangeOuter,
   rangeFull,
   omegaWorkingPrecisionParametric,
   omegaAccuracyGoalParametric,
   omegaPrecisionGoalParametric,
   omegaPrecisionGoalInnerParametric,
   omegaMaxRecursionParametric
   },

  omegaWorkingPrecisionParametric = 30;
  omegaAccuracyGoalParametric = 12;
  omegaPrecisionGoalParametric = 12;
  omegaPrecisionGoalInnerParametric = 4; (* Analytical strategy: relaxed tolerance in the flat region *)
  omegaMaxRecursionParametric = 30;

  ClearAll[
   omegaIntegrandAParametric,
   omegaIntegrandBParametric
   ];

  omegaIntegrandAParametric[R_?NumericQ] :=
   f0Obj["F0"][R]*Sinh[pb["Psi"][R]]*R;

  omegaIntegrandBParametric[R_?NumericQ] :=
   f1Obj["F1"][R]*Sinh[pb["Psi"][R]]*R;

  If[
   deltaInput > 1,
   rSplit = 1 - 1/deltaInput;
   rangeInner = {R, 0, rSplit};
   rangeOuter = {R, rSplit, 1};

   (* Explicit partitioned evaluations *)
   AParametric =
    NIntegrate[omegaIntegrandAParametric[R], Evaluate[rangeInner],
     WorkingPrecision -> omegaWorkingPrecisionParametric, AccuracyGoal -> omegaAccuracyGoalParametric,
     PrecisionGoal -> omegaPrecisionGoalInnerParametric, MaxRecursion -> omegaMaxRecursionParametric,
     Method -> {"GlobalAdaptive", "SymbolicProcessing" -> 0}] +
    NIntegrate[omegaIntegrandAParametric[R], Evaluate[rangeOuter],
     WorkingPrecision -> omegaWorkingPrecisionParametric, AccuracyGoal -> omegaAccuracyGoalParametric,
     PrecisionGoal -> omegaPrecisionGoalParametric, MaxRecursion -> omegaMaxRecursionParametric,
     Method -> {"GlobalAdaptive", "SymbolicProcessing" -> 0}];

   BParametric =
    NIntegrate[omegaIntegrandBParametric[R], Evaluate[rangeInner],
     WorkingPrecision -> omegaWorkingPrecisionParametric, AccuracyGoal -> omegaAccuracyGoalParametric,
     PrecisionGoal -> omegaPrecisionGoalInnerParametric, MaxRecursion -> omegaMaxRecursionParametric,
     Method -> {"GlobalAdaptive", "SymbolicProcessing" -> 0}] +
    NIntegrate[omegaIntegrandBParametric[R], Evaluate[rangeOuter],
     WorkingPrecision -> omegaWorkingPrecisionParametric, AccuracyGoal -> omegaAccuracyGoalParametric,
     PrecisionGoal -> omegaPrecisionGoalParametric, MaxRecursion -> omegaMaxRecursionParametric,
     Method -> {"GlobalAdaptive", "SymbolicProcessing" -> 0}];
   ,
   (* If delta <= 1, the curve is smooth over the entire domain *)
   rSplit = None;
   rangeFull = {R, 0, 1};

   AParametric =
    NIntegrate[omegaIntegrandAParametric[R], Evaluate[rangeFull],
     WorkingPrecision -> omegaWorkingPrecisionParametric, AccuracyGoal -> omegaAccuracyGoalParametric,
     PrecisionGoal -> omegaPrecisionGoalParametric, MaxRecursion -> omegaMaxRecursionParametric,
     Method -> {"GlobalAdaptive", "SymbolicProcessing" -> 0}];

   BParametric =
    NIntegrate[omegaIntegrandBParametric[R], Evaluate[rangeFull],
     WorkingPrecision -> omegaWorkingPrecisionParametric, AccuracyGoal -> omegaAccuracyGoalParametric,
     PrecisionGoal -> omegaPrecisionGoalParametric, MaxRecursion -> omegaMaxRecursionParametric,
     Method -> {"GlobalAdaptive", "SymbolicProcessing" -> 0}];
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

  OmegaParametric = AParametric/denominatorParametric;

  <|
   "A" -> AParametric,
   "B" -> BParametric,
   "Denominator" -> denominatorParametric,
   "Omega" -> OmegaParametric,
   "Delta" -> deltaInput,
   "SplitPoint" -> rSplit
   |>
  ];
