(* ============================================================ *)
(* Phi_av_parametric.wl                                         *)
(* Parametric implementation of Phi_av                          *)
(* ============================================================ *)
(*
   Computes

       Phi_av =
         -(1/2) Lambda delta^2 G Omega

   from previously validated parametric quantities.

   No new differential equation or numerical integration
   is introduced here.
*)
(* ============================================================ *)


ClearAll[ComputePhiAvParametric];


ComputePhiAvParametric[
   lambdaInput_?NumericQ,
   deltaInput_?NumericQ,
   gObj_Association,
   omegaObj_Association
   ] :=
 Module[
  {
   LambdaLoc,
   deltaLoc,
   GLoc,
   OmegaLoc,
   dPhiDZParametric,
   PhiAvParametric
   },


  (* ========================================================= *)
  (* INPUTS FROM THE VALIDATED PARAMETRIC CHAIN                *)
  (* ========================================================= *)

  LambdaLoc =
   lambdaInput;

  deltaLoc =
   deltaInput;

  GLoc =
   gObj["G"];

  OmegaLoc =
   omegaObj["Omega"];


  (* ========================================================= *)
  (* STREAMING-POTENTIAL GRADIENT                              *)
  (* ========================================================= *)

  dPhiDZParametric =
   -LambdaLoc*
    deltaLoc^2*
    GLoc*
    OmegaLoc;


  (* ========================================================= *)
  (* AXIAL AVERAGE                                             *)
  (*                                                         *)
  (* Phi(Z) = (dPhi/dZ) Z, with Phi(0)=0                     *)
  (*                                                         *)
  (* Phi_av = Integral_0^1 Phi(Z) dZ                          *)
  (*        = (1/2) dPhi/dZ                                   *)
  (* ========================================================= *)

  PhiAvParametric =
   dPhiDZParametric/2;


  (* ========================================================= *)
  (* RETURN OBJECT                                             *)
  (* ========================================================= *)

  <|
   "Lambda" -> LambdaLoc,
   "Delta" -> deltaLoc,
   "G" -> GLoc,
   "Omega" -> OmegaLoc,
   "dPhiDZ" -> dPhiDZParametric,
   "PhiAv" -> PhiAvParametric
   |>
  ];