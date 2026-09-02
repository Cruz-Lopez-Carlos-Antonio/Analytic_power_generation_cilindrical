(* ============================================================ *)
(* Code developed by Sánchez Lozano, G., Cruz-López C.-A., and  *)
(* F. Méndez, for the research:                                 *)
(* "Power generation in cylindrical microchannels with high     *)
(* surface zeta potential: a survey on viscoelectric effects    *)
(* over streaming potentials and efficiency"                    *)
(* Year: 2026                                                   *)
(* ------------------------------------------------------------ *)
(* SemianalyticalParametricSolver.wl                            *)
(* General semianalytical multiparametric solver                *)
(* ============================================================ *)
(*
   PURPOSE
   -------

   General wrapper around the already validated parametric chain.

   It does NOT rewrite the numerical logic of the physical modules.
   It only resolves the requested input parameters and propagates
   them through the existing validated chain.

   Usage examples
   --------------

   1) Base values from Parameters.wl:

      sol = SemianalyticalParametricSolver[delta];

   2) External PsiS:

      sol = SemianalyticalParametricSolver[
         16,
         <|"PsiS" -> -13/10|>
      ];

   3) External omega:

      sol = SemianalyticalParametricSolver[
         16,
         <|"omega" -> 2/1000|>
      ];

   4) Forced Lambda sensitivity study:

      sol = SemianalyticalParametricSolver[
         16,
         <|"Lambda" -> 45/100|>
      ];

   5) Forced PiD sensitivity study:

      sol = SemianalyticalParametricSolver[
         16,
         <|"PiD" -> 150|>
      ];

   AUTOMATIC SEMANTICS
   -------------------

      "PsiS"   -> Automatic   uses PsiS from Parameters.wl
      "omega"  -> Automatic   uses omega from Parameters.wl
      "Lambda" -> Automatic   computes ComputeLambda[Fcc]
      "PiD"    -> Automatic   uses alpha xi^2

   EXACT EXTERNAL INPUT RULE
   -------------------------

   Any EXPLICIT value supplied through the configuration Association
   must be an exact number (integer, rational, etc.).

   Examples:
      "PsiS"   -> -13/10
      "omega"  -> 2/1000
      "Lambda" -> 45/100
      "PiD"    -> 150

   Do not supply machine-precision decimals such as -1.30, 0.002,
   0.45 or 150.  The solver deliberately does NOT repair such inputs
   with SetPrecision, because doing so would assign artificial
   precision to information that was never present.

   IMPORTANT NUMERICAL PRINCIPLES
   ------------------------------

   - No artificial SetPrecision is introduced here.
   - Existing WorkingPrecision / AccuracyGoal / PrecisionGoal /
     integration methods remain inside their validated modules.
   - Omega is computed by SolveOmegaParametric[..., deltaInput],
     whose validated implementation contains the adaptive split

         Rsplit = 1 - 1/delta

     for delta > 1.

   - This core solver deliberately does NOT wrap the complete
     evaluation in Check[..., $Failed].
   - It also deliberately does NOT apply Quiet globally.

     Therefore numerical messages remain available during
     individual diagnostic runs.

     Robust classification of large table evaluations must be
     performed externally by inspecting the returned Association
     and the requested numerical observable.

   CURRENT OUTPUT
   --------------

   The present version returns PhiAv and the useful intermediate
   objects. Eta will be incorporated after its definitive module
   is supplied and validated.
*)
(* ============================================================ *)


(* ============================================================ *)
(* 1. MODULE DIRECTORY                                          *)
(* ============================================================ *)

semianalyticalModuleDirectory =
   DirectoryName[$InputFileName];

If[
   semianalyticalModuleDirectory === "",
   Print[
      "ERROR: Run this file with Get[\"full_path\\SemianalyticalParametricSolver.wl\"]."
   ];
   Abort[];
];


(* ============================================================ *)
(* 2. LOAD BASE PARAMETERS                                      *)
(* ============================================================ *)

Get[
   FileNameJoin[
      {
         semianalyticalModuleDirectory,
         "Parameters.wl"
      }
   ]
];


(* ============================================================ *)
(* 3. LOAD VALIDATED PARAMETRIC MODULES                         *)
(* ============================================================ *)

Scan[
   Get[
      FileNameJoin[
         {
            semianalyticalModuleDirectory,
            #
         }
      ]
   ] &,

   {
      "PoissonBoltzmann_parametric.wl",
      "F_cc_parametric.wl",
      "Lambda_parameter_parametric.wl",
      "MR_parametric.wl",
      "HR_parametric.wl",
      "F0_parametric.wl",
      "F1_parametric.wl",
      "Omega_parameter_parametric.wl",
      "F_parametric.wl",
      "G_parametric.wl"
   }
];


(* ============================================================ *)
(* 4. MESSAGES                                                  *)
(* ============================================================ *)

ClearAll[
   SemianalyticalParametricSolver,
   exactExternalNumberQ
];

SemianalyticalParametricSolver::arg =
   "The second argument must be an Association.";

SemianalyticalParametricSolver::key =
   "Unknown configuration key(s): `1`. Allowed keys are \"PsiS\", \"omega\", \"Lambda\", and \"PiD\".";

SemianalyticalParametricSolver::num =
   "Resolved value for \"`1`\" is not numeric: `2`.";

SemianalyticalParametricSolver::inexact =
   "The explicit value supplied for \"`1`\" is inexact (`2`). Use an exact integer or rational value instead (for example, 45/100 instead of 0.45). The solver will not apply SetPrecision to an inexact external input.";

SemianalyticalParametricSolver::delta =
   "deltaInput must be a positive numeric value. Received: `1`.";


(* ============================================================ *)
(* 5. EXACT-INPUT TEST                                          *)
(* ============================================================ *)

exactExternalNumberQ[x_] :=
   NumericQ[x] &&
   Precision[x] === Infinity;


(* ============================================================ *)
(* 6. ONE-ARGUMENT CONVENIENCE FORM                             *)
(* ============================================================ *)

SemianalyticalParametricSolver[
   deltaInput_?NumericQ
   ] :=
   SemianalyticalParametricSolver[
      deltaInput,
      <||>
   ];


(* ============================================================ *)
(* 7. GENERAL MULTIPARAMETRIC SOLVER                            *)
(* ============================================================ *)

SemianalyticalParametricSolver[
   deltaInput_?NumericQ,
   config_Association
   ] :=
 Module[
   {
      allowedKeys,
      unknownKeys,
      defaults,
      cfg,

      psiSRequested,
      omegaRequested,
      lambdaRequested,
      piDRequested,

      psiSValue,
      omegaValue,
      lambdaValue,
      piDValue,

      psiSMode,
      omegaMode,
      lambdaMode,
      piDMode,

      pbObj,
      hObj,
      fccValue,
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


   (* --------------------------------------------------------- *)
   (* INPUT DOMAIN                                              *)
   (* --------------------------------------------------------- *)

   If[
      deltaInput <= 0,
      Message[
         SemianalyticalParametricSolver::delta,
         deltaInput
      ];
      Return[$Failed]
   ];


   (* --------------------------------------------------------- *)
   (* CONFIGURATION                                             *)
   (* --------------------------------------------------------- *)

   allowedKeys = {
      "PsiS",
      "omega",
      "Lambda",
      "PiD"
   };

   unknownKeys =
      Complement[
         Keys[config],
         allowedKeys
      ];

   If[
      unknownKeys =!= {},
      Message[
         SemianalyticalParametricSolver::key,
         unknownKeys
      ];
      Return[$Failed]
   ];


   defaults =
      <|
         "PsiS" -> Automatic,
         "omega" -> Automatic,
         "Lambda" -> Automatic,
         "PiD" -> Automatic
      |>;

   cfg =
      Join[
         defaults,
         config
      ];


   psiSRequested   = cfg["PsiS"];
   omegaRequested  = cfg["omega"];
   lambdaRequested = cfg["Lambda"];
   piDRequested    = cfg["PiD"];


   (* --------------------------------------------------------- *)
   (* VALIDATE EXPLICIT EXTERNAL INPUTS                         *)
   (* --------------------------------------------------------- *)

   If[
      psiSRequested =!= Automatic &&
      ! exactExternalNumberQ[psiSRequested],
      Message[
         SemianalyticalParametricSolver::inexact,
         "PsiS",
         psiSRequested
      ];
      Return[$Failed]
   ];

   If[
      omegaRequested =!= Automatic &&
      ! exactExternalNumberQ[omegaRequested],
      Message[
         SemianalyticalParametricSolver::inexact,
         "omega",
         omegaRequested
      ];
      Return[$Failed]
   ];

   If[
      lambdaRequested =!= Automatic &&
      ! exactExternalNumberQ[lambdaRequested],
      Message[
         SemianalyticalParametricSolver::inexact,
         "Lambda",
         lambdaRequested
      ];
      Return[$Failed]
   ];

   If[
      piDRequested =!= Automatic &&
      ! exactExternalNumberQ[piDRequested],
      Message[
         SemianalyticalParametricSolver::inexact,
         "PiD",
         piDRequested
      ];
      Return[$Failed]
   ];


   (* --------------------------------------------------------- *)
   (* RESOLVE PsiS                                             *)
   (* --------------------------------------------------------- *)

   If[
      psiSRequested === Automatic,

      psiSValue = PsiS;
      psiSMode = "Automatic",

      psiSValue = psiSRequested;
      psiSMode = "Explicit"
   ];

   If[
      ! NumericQ[psiSValue],
      Message[
         SemianalyticalParametricSolver::num,
         "PsiS",
         psiSValue
      ];
      Return[$Failed]
   ];


   (* --------------------------------------------------------- *)
   (* RESOLVE omega                                            *)
   (* --------------------------------------------------------- *)

   If[
      omegaRequested === Automatic,

      omegaValue = omega;
      omegaMode = "Automatic",

      omegaValue = omegaRequested;
      omegaMode = "Explicit"
   ];

   If[
      ! NumericQ[omegaValue],
      Message[
         SemianalyticalParametricSolver::num,
         "omega",
         omegaValue
      ];
      Return[$Failed]
   ];


   (* --------------------------------------------------------- *)
   (* Poisson-Boltzmann                                        *)
   (* --------------------------------------------------------- *)

   pbObj =
      SolvePB[
         deltaInput,
         psiSValue
      ];


   (* --------------------------------------------------------- *)
   (* H(R)                                                      *)
   (* --------------------------------------------------------- *)

   hObj =
      SolveHParametric[
         pbObj
      ];


   (* --------------------------------------------------------- *)
   (* F_cc                                                      *)
   (* --------------------------------------------------------- *)

   fccValue =
      ComputeFcc[
         pbObj
      ];


   (* --------------------------------------------------------- *)
   (* RESOLVE Lambda                                           *)
   (* --------------------------------------------------------- *)

   If[
      lambdaRequested === Automatic,

      lambdaValue =
         ComputeLambda[
            fccValue
         ];
      lambdaMode = "Automatic",

      lambdaValue = lambdaRequested;
      lambdaMode = "Explicit"
   ];

   If[
      ! NumericQ[lambdaValue],
      Message[
         SemianalyticalParametricSolver::num,
         "Lambda",
         lambdaValue
      ];
      Return[$Failed]
   ];


   (* --------------------------------------------------------- *)
   (* RESOLVE PiD                                              *)
   (* --------------------------------------------------------- *)

   If[
      piDRequested === Automatic,

      piDValue =
         alpha xi^2;
      piDMode = "Automatic",

      piDValue = piDRequested;
      piDMode = "Explicit"
   ];

   If[
      ! NumericQ[piDValue],
      Message[
         SemianalyticalParametricSolver::num,
         "PiD",
         piDValue
      ];
      Return[$Failed]
   ];


   (* --------------------------------------------------------- *)
   (* F0(R)                                                     *)
   (* --------------------------------------------------------- *)

   f0Obj =
      SolveF0Parametric[
         pbObj,
         omegaValue,
         piDValue
      ];


   (* --------------------------------------------------------- *)
   (* F1(R)                                                     *)
   (* --------------------------------------------------------- *)

   f1Obj =
      SolveF1Parametric[
         pbObj,
         hObj,
         omegaValue,
         lambdaValue,
         deltaInput
      ];


   (* --------------------------------------------------------- *)
   (* Omega                                                     *)
   (* IMPORTANT: deltaInput is the fourth argument.             *)
   (* The validated Omega module applies Rsplit for delta > 1.  *)
   (* --------------------------------------------------------- *)

   omegaObj =
      SolveOmegaParametric[
         pbObj,
         f0Obj,
         f1Obj,
         deltaInput
      ];

   OmegaValue =
      omegaObj["Omega"];


   (* --------------------------------------------------------- *)
   (* F(R)                                                      *)
   (* --------------------------------------------------------- *)

   fObj =
      BuildFParametric[
         f0Obj,
         f1Obj,
         omegaObj
      ];


   (* --------------------------------------------------------- *)
   (* G                                                         *)
   (* --------------------------------------------------------- *)

   gObj =
      ComputeGParametric[
         fObj
      ];

   GValue =
      gObj["G"];


   (* --------------------------------------------------------- *)
   (* dPhi/dZ                                                   *)
   (* --------------------------------------------------------- *)

   dPhiDZValue =
      -lambdaValue*
       deltaInput^2*
       GValue*
       OmegaValue;


   (* --------------------------------------------------------- *)
   (* Phi_av                                                    *)
   (* --------------------------------------------------------- *)

   PhiAvValue =
      dPhiDZValue/2;


   (* --------------------------------------------------------- *)
   (* RETURN COMPLETE SOLUTION OBJECT                           *)
   (* --------------------------------------------------------- *)

   <|
      "Delta" -> deltaInput,

      "PsiS" -> psiSValue,
      "omega" -> omegaValue,
      "Lambda" -> lambdaValue,
      "PiD" -> piDValue,

      "PsiSMode" -> psiSMode,
      "omegaMode" -> omegaMode,
      "LambdaMode" -> lambdaMode,
      "PiDMode" -> piDMode,

      "Fcc" -> fccValue,

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


(* ============================================================ *)
(* 8. FALLBACK FOR NON-ASSOCIATION SECOND ARGUMENT              *)
(* ============================================================ *)

SemianalyticalParametricSolver[
   deltaInput_?NumericQ,
   other_
   ] :=
 (
   Message[
      SemianalyticalParametricSolver::arg
   ];
   $Failed
 );
