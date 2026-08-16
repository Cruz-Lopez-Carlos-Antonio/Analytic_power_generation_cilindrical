(* ============================================================ *)
(* eta_verification.wl                                         *)
(* Verification of the electrokinetic efficiency eta           *)
(* ============================================================ *)

moduleDirectory = DirectoryName[$InputFileName];

Get[FileNameJoin[{moduleDirectory, "Lambda_parameter.wl"}]];
Get[FileNameJoin[{moduleDirectory, "Omega_parameter.wl"}]];
Get[FileNameJoin[{moduleDirectory, "G.wl"}]];


(* ============================================================ *)
(* CHARACTERISTIC POTENTIAL SCALE                              *)
(* ============================================================ *)

phiC =
  (mu Jw l)/
  (epsilon zeta delta^2);


(* ============================================================ *)
(* CHARACTERISTIC EFFICIENCY eta_c                             *)
(* Equation (63)                                               *)
(* ============================================================ *)

etaC1 =
  ((1/4) sigmaInf Fcc (phiC^2/l^2))/
  (Jw (PosDb/l));


(* Algebraically simplified form                               *)

etaC2 =
  (sigmaInf Fcc mu^2 Jw l)/
  (4 PosDb epsilon^2 zeta^2 delta^4);


(* ============================================================ *)
(* STREAMING-POTENTIAL GRADIENT                                *)
(* ============================================================ *)

dPhiDZ =
  -Lambda delta^2 G Omega;


(* ============================================================ *)
(* EFFICIENCY                                                  *)
(* Equation (62)                                               *)
(* ============================================================ *)

eta =
  etaC1 (dPhiDZ^2)/(-G);


(* ============================================================ *)
(* VERIFICATION OUTPUT                                         *)
(* ============================================================ *)

Print["=============================================="];
Print["EFFICIENCY VERIFICATION"];
Print["=============================================="];

Print["delta       = ", N[delta, 18]];
Print["Lambda      = ", N[Lambda, 18]];
Print["Omega       = ", N[Omega, 18]];
Print["G           = ", N[G, 18]];

Print[""];
Print["sigmaInf    = ", N[sigmaInf, 18]];
Print["Fcc         = ", N[Fcc, 18]];
Print["phiC        = ", N[phiC, 18]];

Print[""];
Print["etaC (Eq.63)       = ", ScientificForm[N[etaC1, 18]]];
Print["etaC (simplified)  = ", ScientificForm[N[etaC2, 18]]];

Print[
  "|etaC1-etaC2|       = ",
  ScientificForm[N[Abs[etaC1 - etaC2], 10]]
];

Print[""];
Print["dPhi/dZ            = ", N[dPhiDZ, 18]];
Print["eta                = ", ScientificForm[N[eta, 18]]];

Print["=============================================="];