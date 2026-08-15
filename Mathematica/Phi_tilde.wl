(* ============================================================ *)
(* PiTilde.wl                                                   *)
(* Dimensionless modified pressure                              *)
(* ============================================================ *)

Module[{localDir = DirectoryName[$InputFileName]},

  Get[FileNameJoin[{localDir, "G.wl"}]];
];

(* ============================================================ *)
(* MODIFIED PRESSURE                                            *)
(* ============================================================ *)

ClearAll[PiTilde];

PiTilde[Z_?NumericQ] := 1 + G*Z;