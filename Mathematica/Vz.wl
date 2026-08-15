(* ============================================================ *)
(* Vz.wl                                                       *)
(* Dimensionless axial velocity V_Z(R)                         *)
(* ============================================================ *)
(*                                                            *)
(* Operational module for                                     *)
(*                                                            *)
(*     V_Z(R) = G F(R)                                        *)
(*                                                            *)
(* and its radial derivative                                  *)
(*                                                            *)
(*     V_Z'(R) = G F'(R).                                     *)
(*                                                            *)
(* This module is intended for downstream use.                 *)
(* It remains silent when imported.                            *)
(*                                                            *)
(* ============================================================ *)

Module[{localDir = DirectoryName[$InputFileName]},
  Get[FileNameJoin[{localDir, "G.wl"}]];
];

(* ============================================================ *)
(* AXIAL VELOCITY                                               *)
(* ============================================================ *)

ClearAll[Vz];

Vz[R_?NumericQ] := G F[R];

(* ============================================================ *)
(* RADIAL DERIVATIVE                                            *)
(* ============================================================ *)

ClearAll[VzPrime];

VzPrime[R_?NumericQ] := G FPrime[R];