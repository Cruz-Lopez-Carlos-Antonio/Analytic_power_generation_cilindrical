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
  omegaPrecisionGoalInnerParametric = 4; (* <-- Magia analítica: tolerancia relajada en zona plana *)
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

   (* Evaluaciones explícitas particionadas *)
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
   (* Si delta <= 1, la curva es suave en todo el dominio *)
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
