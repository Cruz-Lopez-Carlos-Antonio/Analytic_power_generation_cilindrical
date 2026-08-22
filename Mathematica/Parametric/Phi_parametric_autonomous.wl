(* ============================================================ *)
(* Phi_av_parametric.wl                                         *)
(* Autonomous parametric solver for Phi_av                      *)
(* ============================================================ *)
(*
   Computes the complete semianalytical parametric chain

       (delta, PsiS)
          -> Psi(R)
          -> H(R)
          -> Fcc
          -> Lambda
          -> M(R)
          -> F0(R), F1(R)
          -> Omega
          -> F(R)
          -> G
          -> dPhi/dZ
          -> Phi_av

   The validated baseline value of omega is taken from
   Parameters.wl.

   Usage:

       sol = PhiAvParametric[deltaValue, psiSValue];

       sol["PhiAv"]

   Other quantities can also be inspected:

       sol["Lambda"]
       sol["Omega"]
       sol["G"]
       sol["dPhiDZ"]

   IMPORTANT:
   No additional SetPrecision manipulations are introduced.
*)
(* ============================================================ *)


(* ============================================================ *)
(* 1. MODULE DIRECTORY                                         *)
(* ============================================================ *)

phiAvModuleDirectory =
  DirectoryName[$InputFileName];


(* ============================================================ *)
(* 2. LOAD BASE PARAMETERS                                     *)
(* ============================================================ *)

Get[
  FileNameJoin[
    {
      phiAvModuleDirectory,
      "Parameters.wl"
    }
  ]
];


(* ============================================================ *)
(* 3. LOAD VALIDATED PARAMETRIC MODULES                         *)
(* ============================================================ *)

Get[
  FileNameJoin[
    {
      phiAvModuleDirectory,
      "PoissonBoltzmann_parametric.wl"
    }
  ]
];


Get[
  FileNameJoin[
    {
      phiAvModuleDirectory,
      "F_cc_parametric.wl"
    }
  ]
];


Get[
  FileNameJoin[
    {
      phiAvModuleDirectory,
      "Lambda_parameter_parametric.wl"
    }
  ]
];


Get[
  FileNameJoin[
    {
      phiAvModuleDirectory,
      "MR_parametric.wl"
    }
  ]
];


Get[
  FileNameJoin[
    {
      phiAvModuleDirectory,
      "HR_parametric.wl"
    }
  ]
];


Get[
  FileNameJoin[
    {
      phiAvModuleDirectory,
      "F0_parametric.wl"
    }
  ]
];


Get[
  FileNameJoin[
    {
      phiAvModuleDirectory,
      "F1_parametric.wl"
    }
  ]
];


Get[
  FileNameJoin[
    {
      phiAvModuleDirectory,
      "Omega_parameter_parametric.wl"
    }
  ]
];


Get[
  FileNameJoin[
    {
      phiAvModuleDirectory,
      "F_parametric.wl"
    }
  ]
];


Get[
  FileNameJoin[
    {
      phiAvModuleDirectory,
      "G_parametric.wl"
    }
  ]
];


(* ============================================================ *)
(* 4. AUTONOMOUS Phi_av SOLVER                                 *)
(* ============================================================ *)

ClearAll[PhiAvParametric];


PhiAvParametric[
   deltaInput_?NumericQ,
   psiSInput_?NumericQ
   ] :=
 Module[
  {
   pbObj,
   hObj,
   fccValue,
   lambdaValue,
   piDValue,
   f0Obj,
   f1Obj,
   omegaObj,
   fObj,
   gObj,
   OmegaValue,
   GValue,
   dPhiDZValue,
   PhiAvValue
   },


  (* ========================================================= *)
  (* POISSON-BOLTZMANN                                        *)
  (* ========================================================= *)

  pbObj =
    SolvePB[
      deltaInput,
      psiSInput
    ];


  (* ========================================================= *)
  (* H(R)                                                      *)
  (* ========================================================= *)

  hObj =
    SolveHParametric[
      pbObj
    ];


  (* ========================================================= *)
  (* F_cc                                                      *)
  (* ========================================================= *)

  fccValue =
    ComputeFcc[
      pbObj
    ];


  (* ========================================================= *)
  (* Lambda                                                    *)
  (* ========================================================= *)

  lambdaValue =
    ComputeLambda[
      fccValue
    ];


  (* ========================================================= *)
  (* Pi_D                                                      *)
  (* ========================================================= *)
  (*
     For the current parametric study, the physical parameters
     entering alpha and xi remain those of Parameters.wl.
  *)

  piDValue =
    alpha xi^2;


  (* ========================================================= *)
  (* F0(R)                                                     *)
  (* ========================================================= *)

  f0Obj =
    SolveF0Parametric[
      pbObj,
      omega,
      piDValue
    ];


  (* ========================================================= *)
  (* F1(R)                                                     *)
  (* ========================================================= *)

  f1Obj =
    SolveF1Parametric[
      pbObj,
      hObj,
      omega,
      lambdaValue,
      deltaInput
    ];


  (* ========================================================= *)
  (* Omega                                                     *)
  (* ========================================================= *)

  omegaObj =
    SolveOmegaParametric[
      pbObj,
      f0Obj,
      f1Obj
    ];


  OmegaValue =
    omegaObj["Omega"];


  (* ========================================================= *)
  (* F(R)                                                      *)
  (* ========================================================= *)

  fObj =
    BuildFParametric[
      f0Obj,
      f1Obj,
      omegaObj
    ];


  (* ========================================================= *)
  (* G                                                         *)
  (* ========================================================= *)

  gObj =
    ComputeGParametric[
      fObj
    ];


  GValue =
    gObj["G"];


  (* ========================================================= *)
  (* dPhi/dZ                                                   *)
  (* ========================================================= *)

  dPhiDZValue =
    -lambdaValue*
     deltaInput^2*
     GValue*
     OmegaValue;


  (* ========================================================= *)
  (* Phi_av                                                    *)
  (* ========================================================= *)

  PhiAvValue =
    dPhiDZValue/2;


  (* ========================================================= *)
  (* RETURN COMPLETE SOLUTION OBJECT                           *)
  (* ========================================================= *)

  <|
    "Delta" -> deltaInput,
    "PsiS" -> psiSInput,

    "Fcc" -> fccValue,
    "Lambda" -> lambdaValue,
    "PiD" -> piDValue,

    "Omega" -> OmegaValue,
    "G" -> GValue,

    "dPhiDZ" -> dPhiDZValue,
    "PhiAv" -> PhiAvValue,

    "PBObject" -> pbObj,
    "HObject" -> hObj,
    "F0Object" -> f0Obj,
    "F1Object" -> f1Obj,
    "OmegaObject" -> omegaObj,
    "FObject" -> fObj,
    "GObject" -> gObj
  |>
];