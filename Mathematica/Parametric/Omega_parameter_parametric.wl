(* ============================================================ *)
(* Omega_parameter_parametric.wl                                *)
(* Parametric implementation of A, B and Omega                  *)
(* ============================================================ *)
(*
   Mirrors the validated Omega_parameter.wl module.

   Definitions:

       A =
         Integral_0^1
         F0(R) Sinh[Psi(R)] R dR

       B =
         Integral_0^1
         F1(R) Sinh[Psi(R)] R dR

       Omega =
         A/(1-B)

   The numerical settings are intentionally kept identical
   to the original operational module.

   The only substantive changes are that Psi, F0 and F1
   are supplied through previously constructed parametric
   solution objects.
*)
(* ============================================================ *)


ClearAll[SolveOmegaParametric];


SolveOmegaParametric[
   pb_Association,
   f0Obj_Association,
   f1Obj_Association
   ] :=
 Module[
  {
   omegaIntegrandAParametric,
   omegaIntegrandBParametric,
   AParametric,
   BParametric,
   denominatorParametric,
   OmegaParametric
   },


  (* ========================================================= *)
  (* NUMERICAL SETTINGS                                        *)
  (* IDENTICAL TO Omega_parameter.wl                            *)
  (* ========================================================= *)

  omegaWorkingPrecisionParametric = 30;
  omegaAccuracyGoalParametric = 12;
  omegaPrecisionGoalParametric = 12;
  omegaMaxRecursionParametric = 30;


  (* ========================================================= *)
  (* INTEGRANDS                                                *)
  (* ========================================================= *)

  ClearAll[
   omegaIntegrandAParametric,
   omegaIntegrandBParametric
   ];


  omegaIntegrandAParametric[R_?NumericQ] :=
   f0Obj["F0"][R]*
    Sinh[
     pb["Psi"][R]
     ]*
    R;


  omegaIntegrandBParametric[R_?NumericQ] :=
   f1Obj["F1"][R]*
    Sinh[
     pb["Psi"][R]
     ]*
    R;


  (* ========================================================= *)
  (* A                                                        *)
  (* ========================================================= *)

  AParametric =
   NIntegrate[
    omegaIntegrandAParametric[R],

    {R, 0, 1},

    WorkingPrecision ->
     omegaWorkingPrecisionParametric,

    AccuracyGoal ->
     omegaAccuracyGoalParametric,

    PrecisionGoal ->
     omegaPrecisionGoalParametric,

    MaxRecursion ->
     omegaMaxRecursionParametric,

    Method -> {
      "GlobalAdaptive",
      "SymbolicProcessing" -> 0
      }
    ];


  (* ========================================================= *)
  (* B                                                        *)
  (* ========================================================= *)

  BParametric =
   NIntegrate[
    omegaIntegrandBParametric[R],

    {R, 0, 1},

    WorkingPrecision ->
     omegaWorkingPrecisionParametric,

    AccuracyGoal ->
     omegaAccuracyGoalParametric,

    PrecisionGoal ->
     omegaPrecisionGoalParametric,

    MaxRecursion ->
     omegaMaxRecursionParametric,

    Method -> {
      "GlobalAdaptive",
      "SymbolicProcessing" -> 0
      }
    ];


  (* ========================================================= *)
  (* OMEGA                                                    *)
  (* ========================================================= *)

  denominatorParametric =
   1 - BParametric;


  If[
   Abs[
     denominatorParametric
     ] < 10^-8,

   Print[
    "WARNING in Omega_parameter_parametric.wl: ",
    "1 - B is numerically small: ",
    N[
     denominatorParametric,
     18
     ]
    ]
   ];


  OmegaParametric =
   AParametric/
    denominatorParametric;


  (* ========================================================= *)
  (* RETURN OBJECT                                            *)
  (* ========================================================= *)

  <|
   "A" -> AParametric,
   "B" -> BParametric,
   "Denominator" -> denominatorParametric,
   "Omega" -> OmegaParametric
   |>
  ];