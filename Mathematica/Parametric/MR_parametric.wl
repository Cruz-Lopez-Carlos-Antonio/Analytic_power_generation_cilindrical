(* ============================================================ *)
(* MR_parametric.wl                                             *)
(* Parametric implementation of M(R)                            *)
(* ============================================================ *)
(*
   Computes

       M(R) = Exp[omega (Psi'(R))^2]

   from a previously constructed parametric Poisson-Boltzmann
   solution object:

       pb = SolvePB[deltaValue, psiSValue];

   The numerical structure intentionally mirrors the original
   MR.wl module as closely as possible.

   IMPORTANT:
   No artificial SetPrecision is applied to omegaInput or RInput.
   Their native precision is preserved.
*)
(* ============================================================ *)


ClearAll[ComputeM];


(* ============================================================ *)
(* PARAMETRIC M(R)                                              *)
(* ============================================================ *)

ComputeM[
   pb_Association,
   omegaInput_?NumericQ,
   RInput_?NumericQ
   ] :=
 Module[
  {
   psiPrimeValue
   },


  (* ---------------------------------------------------------- *)
  (* DOMAIN CHECK                                               *)
  (* ---------------------------------------------------------- *)

  If[
   RInput < 0 || RInput > 1,

   Return[
    Indeterminate
    ]
   ];


  (* ---------------------------------------------------------- *)
  (* Psi'(R) FROM PARAMETRIC PB OBJECT                          *)
  (* ---------------------------------------------------------- *)

  psiPrimeValue =
   pb["PsiPrime"][RInput];


  (* ---------------------------------------------------------- *)
  (* M(R)                                                       *)
  (* ---------------------------------------------------------- *)

  Exp[
   omegaInput*
    psiPrimeValue^2
   ]
  ];