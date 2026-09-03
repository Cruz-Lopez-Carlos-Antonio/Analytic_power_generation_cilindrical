(* ============================================================ *)
(* Code developed by Sánchez Lozano, G., Cruz-López C.-A., and  *)
(* F. Méndez, for the research:                                 *)
(* "Power generation in cylindrical microchannels with high     *)
(* surface zeta potential: a survey on viscoelectric effects    *)
(* over streaming potentials and efficiency"                    *)
(* Year: 2026                                                   *)
(* ------------------------------------------------------------ *)
(* F0_parametric.wl                                             *)
(* Parametric implementation of F0(R)                           *)
(* ============================================================ *)
(*
   This module mirrors the validated F0.wl implementation.

   Original equation:

       F0'(R) = PiD R / (2 M(R)),
       F0(1)  = 0.

   The numerical strategy is intentionally kept identical
   to the original implementation.

   The only substantive change is that M(R) is supplied
   through the already validated parametric chain.
*)
(* ============================================================ *)


ClearAll[SolveF0Parametric];


SolveF0Parametric[
   pb_Association,
   omegaInput_?NumericQ,
   piDInput_?NumericQ
   ] :=
 Module[
  {
   f0,
   R,
   f0Sol,
   F0Parametric,
   F0PrimeParametric,
   piDLoc
   },


  (* ========================================================= *)
  (* LOCAL PARAMETER                                           *)
  (* ========================================================= *)

  piDLoc = piDInput;


  (* ========================================================= *)
  (* AUXILIARY ODE                                             *)
  (* ========================================================= *)

  f0Sol =
    NDSolveValue[
      {
        f0'[R] ==
          piDLoc R/
           (
            2 ComputeM[
              pb,
              omegaInput,
              R
            ]
           ),

        f0[1] == 0
      },

      f0,

      {R, 0, 1},

      WorkingPrecision -> 30,
      AccuracyGoal -> 14,
      PrecisionGoal -> 14,

      Method -> {
        "TimeIntegration" -> {
          "ExplicitRungeKutta",
          "DifferenceOrder" -> 8
        }
      }
    ];


  (* ========================================================= *)
  (* PUBLIC FUNCTION                                            *)
  (* ========================================================= *)

  F0Parametric[Rin_?NumericQ] :=
    Which[

      0 <= Rin <= 1,
      f0Sol[Rin],

      True,
      Indeterminate
    ];


  (* ========================================================= *)
  (* PUBLIC DERIVATIVE — EXACT OPERATIONAL RHS                  *)
  (* ========================================================= *)

  F0PrimeParametric[Rin_?NumericQ] :=
    Which[

      0 <= Rin <= 1,
      piDLoc Rin/
       (
        2 ComputeM[
          pb,
          omegaInput,
          Rin
        ]
       ),

      True,
      Indeterminate
    ];


  (* ========================================================= *)
  (* RETURN OBJECT                                              *)
  (* ========================================================= *)

  <|
    "F0" -> F0Parametric,
    "F0Prime" -> F0PrimeParametric,
    "Solution" -> f0Sol,
    "PiD" -> piDLoc
  |>
];
