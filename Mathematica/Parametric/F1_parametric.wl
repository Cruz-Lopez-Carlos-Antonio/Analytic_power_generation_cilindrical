(* ============================================================ *)
(* Code developed by Sánchez Lozano, G., Cruz-López C.-A., and  *)
(* F. Méndez, for the research:                                 *)
(* "Power generation in cylindrical microchannels with high     *)
(* surface zeta potential: a survey on viscoelectric effects    *)
(* over streaming potentials and efficiency"                    *)
(* Year: 2026                                                   *)
(* ------------------------------------------------------------ *)
(* F1_parametric.wl                                             *)
(* Strict parametric clone of validated F1.wl                   *)
(* ============================================================ *)

ClearAll[SolveF1Parametric];

SolveF1Parametric[
   pb_Association,
   hObj_Association,
   omegaInput_?NumericQ,
   lambdaInput_?NumericQ,
   deltaInput_?NumericQ
   ] :=
 Module[
  {
   HRParam,
   MRParam,
   QF1Param,
   f1,
   f1Sol,
   F1Param,
   F1PrimeParam,
   wpF1 = 30,
   agF1 = 14,
   pgF1 = 14,
   epsQ = 10^-8,
   LambdaLoc,
   deltaLoc
   },

  LambdaLoc = lambdaInput;
  deltaLoc = deltaInput;


  (* ---------------------------------------------------------- *)
  (* Local analogues of original HR and MR                      *)
  (* ---------------------------------------------------------- *)

  HRParam[Rin_?NumericQ] :=
    hObj["H"][Rin];

  MRParam[Rin_?NumericQ] :=
    ComputeM[
      pb,
      omegaInput,
      Rin
    ];


  (* ---------------------------------------------------------- *)
  (* Auxiliary function QF1(R)                                  *)
  (* Structure intentionally identical to original F1.wl        *)
  (* ---------------------------------------------------------- *)

  QF1Param[Rin_?NumericQ] := Which[
      Rin == 0, 0,

      0 < Rin < epsQ,
      (Sinh[pb["CenterPotential"]]/2) Rin,

      epsQ <= Rin <= 1,
      HRParam[Rin]/(Rin MRParam[Rin]),

      True,
      Indeterminate
  ];


  (* ---------------------------------------------------------- *)
  (* Auxiliary ODE                                              *)
  (* ---------------------------------------------------------- *)

  f1Sol = NDSolveValue[
      {
          f1'[R] ==
            LambdaLoc deltaLoc^2 QF1Param[R],

          f1[1] == 0
      },

      f1,

      {R, 0, 1},

      WorkingPrecision -> wpF1,
      AccuracyGoal -> agF1,
      PrecisionGoal -> pgF1,

      Method -> {
          "TimeIntegration" -> {
              "ExplicitRungeKutta",
              "DifferenceOrder" -> 8
          }
      }
  ];


  (* ---------------------------------------------------------- *)
  (* Public function                                            *)
  (* ---------------------------------------------------------- *)

  F1Param[Rin_?NumericQ] := Which[
      0 <= Rin <= 1,
      f1Sol[Rin],

      True,
      Indeterminate
  ];


  (* ---------------------------------------------------------- *)
  (* Public derivative — operational exact RHS                  *)
  (* ---------------------------------------------------------- *)

  F1PrimeParam[Rin_?NumericQ] := Which[
      0 <= Rin <= 1,
      LambdaLoc deltaLoc^2 QF1Param[Rin],

      True,
      Indeterminate
  ];


  <|
    "F1" -> F1Param,
    "F1Prime" -> F1PrimeParam,
    "QF1" -> QF1Param,
    "Solution" -> f1Sol
  |>
];
