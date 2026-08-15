(* ============================================================ *)
(* Pi_RZ.wl                                                     *)
(* Total dimensionless pressure                                 *)
(* ============================================================ *)

ClearAll[Pressure];

(* ============================================================ *)
(* LOAD PREVIOUS VALIDATED MODULE                              *)
(* ============================================================ *)

Get[
  FileNameJoin[
    {DirectoryName[$InputFileName], "Phi_tilde.wl"}
  ]
];

(* ============================================================ *)
(* TOTAL DIMENSIONLESS PRESSURE                                 *)
(* ============================================================ *)
(*                                                              *)
(* Mathematical notation:                                       *)
(*                                                              *)
(*   Pi(R,Z) = PiTilde(Z) + cosh(Psi(R))                        *)
(*                                                              *)
(* PsiPB[R] is already available through the dependency chain   *)
(* used by Phi_tilde -> G -> F -> ...                           *)
(*                                                              *)
(* ============================================================ *)

Pressure[R_?NumericQ, Z_?NumericQ] :=
  PiTilde[Z] + Cosh[PsiPB[R]];