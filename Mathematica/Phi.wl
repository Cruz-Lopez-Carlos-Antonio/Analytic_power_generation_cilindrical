(* ============================================================ *)
(* Phi.wl                                                       *)
(* Dimensionless streaming potential                            *)
(* ============================================================ *)

ClearAll["Global`*"];

Module[{localDir = DirectoryName[$InputFileName]},

  Get[FileNameJoin[{localDir, "G.wl"}]];
  Get[FileNameJoin[{localDir, "Omega_parameter.wl"}]];
  Get[FileNameJoin[{localDir, "Lambda_parameter.wl"}]];
];

(* ------------------------------------------------------------ *)
(* Streaming-potential gradient                                 *)
(* ------------------------------------------------------------ *)

ClearAll[PhiGradient];

PhiGradient = -Lambda*delta^2*G*Omega;

(* ------------------------------------------------------------ *)
(* Dimensionless streaming potential                            *)
(* ------------------------------------------------------------ *)

ClearAll[Phi];

Phi[Z_?NumericQ] := PhiGradient*Z;

(* ------------------------------------------------------------ *)
(* Average streaming potential over 0 <= Z <= 1                 *)
(* ------------------------------------------------------------ *)

ClearAll[PhiAv];

PhiAv = PhiGradient/2;
