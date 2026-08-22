(* ============================================================ *)
(* F_parametric.wl                                              *)
(* Parametric implementation of F(R)                            *)
(* ============================================================ *)
(*
   Reconstructs

       F(R) = F0(R) + Omega F1(R)

   from previously validated parametric objects.

   Omega is NOT supplied as an arbitrary external parameter.
   It is taken from the parametric Omega object obtained from

       Omega = A/(1-B).

   The derivative is reconstructed operationally as

       F'(R) = F0'(R) + Omega F1'(R).
*)
(* ============================================================ *)


ClearAll[BuildFParametric];


BuildFParametric[
   f0Obj_Association,
   f1Obj_Association,
   omegaObj_Association
   ] :=
 Module[
  {
   OmegaLoc,
   FParametric,
   FPrimeParametric
   },


  (* ========================================================= *)
  (* SELF-CONSISTENT OMEGA                                    *)
  (* ========================================================= *)

  OmegaLoc =
   omegaObj["Omega"];


  (* ========================================================= *)
  (* F(R)                                                     *)
  (* ========================================================= *)

  ClearAll[FParametric];

  FParametric[Rin_?NumericQ] :=
   Which[

    0 <= Rin <= 1,

    f0Obj["F0"][Rin]
     +
     OmegaLoc*
      f1Obj["F1"][Rin],

    True,
    Indeterminate
    ];


  (* ========================================================= *)
  (* F'(R) — OPERATIONAL DERIVATIVE                           *)
  (* ========================================================= *)

  ClearAll[FPrimeParametric];

  FPrimeParametric[Rin_?NumericQ] :=
   Which[

    0 <= Rin <= 1,

    f0Obj["F0Prime"][Rin]
     +
     OmegaLoc*
      f1Obj["F1Prime"][Rin],

    True,
    Indeterminate
    ];


  (* ========================================================= *)
  (* RETURN OBJECT                                            *)
  (* ========================================================= *)

  <|
   "F" -> FParametric,
   "FPrime" -> FPrimeParametric,
   "Omega" -> OmegaLoc
   |>
  ];