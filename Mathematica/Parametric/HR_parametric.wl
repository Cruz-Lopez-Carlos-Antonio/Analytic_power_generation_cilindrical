(* ========================================================= *)
(* HR_parametric.wl                                          *)
(* ========================================================= *)
(*
   Parametric version of the validated HR.wl module.

   Original equation:

       H'(R) = R Sinh[Psi(R)],
       H(0)  = 0.

   The numerical strategy is intentionally kept identical
   to the original HR.wl implementation.

   The ONLY substantive change is that Psi(R) is obtained
   from a previously constructed parametric PB object:

       pb = SolvePB[deltaValue, psiSValue];

*)
(* ========================================================= *)


ClearAll[SolveHParametric];


SolveHParametric[pb_Association] :=
 Module[
  {
   h,
   R,
   hSolution,
   HRParametric
   },


  (* ======================================================= *)
  (* AUXILIARY ODE                                           *)
  (* ======================================================= *)

  hSolution =
    NDSolveValue[
      {
        h'[R] ==
          R Sinh[
            pb["Psi"][R]
          ],

        h[0] == 0
      },

      h,

      {
        R,
        0,
        1
      },

      WorkingPrecision -> 30,
      AccuracyGoal -> 18,
      PrecisionGoal -> 18,

      MaxSteps -> Infinity
    ];


  (* ======================================================= *)
  (* PUBLIC FUNCTION                                         *)
  (* ======================================================= *)

  HRParametric[Rin_?NumericQ] :=
    Which[

      Rin == 0,
      0,

      0 < Rin <= 1,
      hSolution[Rin],

      True,
      Indeterminate
    ];


  (* ======================================================= *)
  (* RETURN OBJECT                                           *)
  (* ======================================================= *)

  <|
    "H" -> HRParametric,
    "Solution" -> hSolution
  |>
];