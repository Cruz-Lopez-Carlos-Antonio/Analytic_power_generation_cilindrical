(* ========================================================= *)
(*
   Computation of the dimensionless parameter

                2 epsilon^2 kappa^2 zetaT^2
       Lambda = -----------------------------
                       mu sigmaInf Fcc

   F_cc.wl loads all previous dependencies required for
   this calculation.
*)
(* ========================================================= *)


(* ========================================================= *)
(* 1. LOAD F_cc MODULE                                       *)
(* ========================================================= *)

moduleDirectory =
  DirectoryName[$InputFileName];

Get[
  FileNameJoin[
    {
      moduleDirectory,
      "F_cc.wl"
    }
  ]
];


(* ========================================================= *)
(* 2. COMPUTE Lambda                                         *)
(* ========================================================= *)

Lambda =
  (2 epsilon^2 kappa^2 zetaT^2)/
  (mu sigmaInf Fcc);
