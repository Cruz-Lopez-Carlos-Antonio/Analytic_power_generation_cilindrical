(* ============================================================ *)
(* Lambda_parameter_parametric.wl                               *)
(* Parametric implementation of Lambda                          *)
(* ============================================================ *)
(*
   Computes the dimensionless parameter

                    2 epsilon^2 kappa^2 zetaT^2
       Lambda = ---------------------------------
                    mu sigmaInf Fcc

   The physical quantities

       epsilon
       kappa
       zetaT
       mu
       sigmaInf

   are inherited from Parameters.wl.

   Fcc is supplied externally from the already validated
   parametric F_cc module.

   IMPORTANT:
   No artificial SetPrecision is applied to fccInput.
   Its native precision is preserved.
*)
(* ============================================================ *)


ClearAll[ComputeLambda];


(* ============================================================ *)
(* PARAMETRIC LAMBDA                                            *)
(* ============================================================ *)

ComputeLambda[
   fccInput_?NumericQ
   ] :=
 (
  2
   epsilon^2
   kappa^2
   zetaT^2
  )/
 (
  mu
   sigmaInf
   fccInput
  );