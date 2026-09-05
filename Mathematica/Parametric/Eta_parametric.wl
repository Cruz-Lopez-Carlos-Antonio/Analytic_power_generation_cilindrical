(* ============================================================ *)
(* Code developed by Sánchez Lozano, G., Cruz-López C.-A., and  *)
(* F. Méndez, for the research:                                 *)
(* "Power generation in cylindrical microchannels with high     *)
(* surface zeta potential: a survey on viscoelectric effects    *)
(* over streaming potentials and efficiency"                    *)
(* Year: 2026                                                   *)
(* ------------------------------------------------------------ *)
(* Eta_parametric.wl                                            *)
(* Parametric computation and internal verification of eta      *)
(* ============================================================ *)
(*
   PURPOSE
   -------

   Compute the electrokinetic conversion efficiency eta from the
   Association returned by SemianalyticalParametricSolver.

   This module deliberately DOES NOT recompute any integral or any
   element of the validated semianalytical chain.  It uses only the
   already available quantities

       Fcc, Omega, G, dPhiDZ, Lambda, Delta, PiD.

   It also performs two independent algebraic evaluations of eta_c
   and two equivalent evaluations of eta.

   DEFINITIONS
   -----------

       phi_c = mu Jw l / (epsilon zeta delta^2)

       eta_c = [(1/4) sigmaInf Fcc (phi_c^2/l^2)] /
               [Jw (PosDb/l)]

       eta   = eta_c (dPhi/dZ)^2 / (-G)

   with

       dPhi/dZ = -Lambda delta^2 G Omega.

   After cancelling the explicit delta^4 factors,

       eta = - C Fcc Lambda^2 G Omega^2,

   where

       C = sigmaInf mu^2 Jw l /
           (4 PosDb epsilon^2 zeta^2).

   IMPORTANT
   ---------

   The symbol names sigmaInf, mu, Jw, l, epsilon, zeta and PosDb
   follow the already validated Parameters.wl / eta verification
   code.  Therefore this module is intended to be loaded after the
   general semianalytical solver (which itself loads Parameters.wl).

   PiD = 0 is a removable singular case for eta in the full chain:

       eta -> 0  as PiD -> 0.

   However, the ordinary general solver computes G = 1/(2 IF), and
   IF vanishes at PiD = 0.  Therefore PiD = 0 should be intercepted
   by the future table generator BEFORE calling the general solver.
*)
(* ============================================================ *)


ClearAll[
   ComputeEtaParametric,
   etaRelativeDifference
];


(* ============================================================ *)
(* MESSAGES                                                     *)
(* ============================================================ *)

ComputeEtaParametric::keys =
   "The supplied solution Association is missing required key(s): `1`.";

ComputeEtaParametric::phys =
   "One or more physical parameters required for eta are not numeric. Check Parameters.wl (sigmaInf, mu, Jw, l, epsilon, zeta, PosDb).";

ComputeEtaParametric::pid0 =
   "PiD = 0 is a removable singular case in the full semianalytical chain. For a parameter study, intercept this case before calling the general solver and assign Eta -> 0 analytically.";

ComputeEtaParametric::verify =
   "The two eta evaluations differ by more than the internal verification tolerance. Relative difference = `1`.";


(* ============================================================ *)
(* RELATIVE DIFFERENCE HELPER                                   *)
(* ============================================================ *)

etaRelativeDifference[x_?NumericQ, y_?NumericQ] :=
 Module[{scale},
   scale = Max[Abs[x], Abs[y]];
   If[scale == 0, 0, Abs[x - y]/scale]
 ];


(* ============================================================ *)
(* MAIN MODULE                                                  *)
(* ============================================================ *)

ComputeEtaParametric[
   sol_Association
   ] :=
 Module[
   {
      requiredKeys,
      missingKeys,

      deltaValue,
      lambdaValue,
      piDValue,
      fccValue,
      OmegaValue,
      GValue,
      dPhiDZValue,

      phiCValue,
      etaCFromDefinition,
      etaCSimplified,
      etaCAbsDifference,
      etaCRelDifference,

      etaFromDefinition,
      etaSimplified,
      etaAbsDifference,
      etaRelDifference,

      prefactorC,
      verificationTolerance,
      etaAgreementQ,
      etaCAgreementQ
   },


   (* --------------------------------------------------------- *)
   (* REQUIRED SOLVER OUTPUT                                    *)
   (* --------------------------------------------------------- *)

   requiredKeys = {
      "Delta",
      "Lambda",
      "PiD",
      "Fcc",
      "Omega",
      "G",
      "dPhiDZ"
   };

   missingKeys =
      Select[
         requiredKeys,
         ! KeyExistsQ[sol, #] &
      ];

   If[
      missingKeys =!= {},
      Message[
         ComputeEtaParametric::keys,
         missingKeys
      ];
      Return[$Failed]
   ];


   (* --------------------------------------------------------- *)
   (* READ SOLVER VALUES                                        *)
   (* --------------------------------------------------------- *)

   deltaValue   = sol["Delta"];
   lambdaValue  = sol["Lambda"];
   piDValue     = sol["PiD"];
   fccValue     = sol["Fcc"];
   OmegaValue   = sol["Omega"];
   GValue       = sol["G"];
   dPhiDZValue  = sol["dPhiDZ"];


   (* --------------------------------------------------------- *)
   (* PHYSICAL-PARAMETER CHECK                                  *)
   (* --------------------------------------------------------- *)

   If[
      ! And @@ (NumericQ /@ {
         sigmaInf,
         mu,
         Jw,
         l,
         epsilon,
         zeta,
         PosDb
      }),
      Message[ComputeEtaParametric::phys];
      Return[$Failed]
   ];


   (* --------------------------------------------------------- *)
   (* REMOVABLE PiD = 0 CASE                                    *)
   (* --------------------------------------------------------- *)

   If[
      TrueQ[piDValue == 0],
      Message[ComputeEtaParametric::pid0]
   ];


   (* --------------------------------------------------------- *)
   (* CHARACTERISTIC POTENTIAL SCALE                            *)
   (* --------------------------------------------------------- *)

   phiCValue =
      (mu Jw l)/
      (epsilon zeta deltaValue^2);


   (* --------------------------------------------------------- *)
   (* eta_c : ROUTE 1 -- ORIGINAL DEFINITION                    *)
   (* --------------------------------------------------------- *)

   etaCFromDefinition =
      ((1/4) sigmaInf fccValue (phiCValue^2/l^2))/
      (Jw (PosDb/l));


   (* --------------------------------------------------------- *)
   (* eta_c : ROUTE 2 -- ALGEBRAICALLY SIMPLIFIED              *)
   (* --------------------------------------------------------- *)

   etaCSimplified =
      (sigmaInf fccValue mu^2 Jw l)/
      (4 PosDb epsilon^2 zeta^2 deltaValue^4);


   etaCAbsDifference =
      Abs[etaCFromDefinition - etaCSimplified];

   etaCRelDifference =
      etaRelativeDifference[
         etaCFromDefinition,
         etaCSimplified
      ];


   (* --------------------------------------------------------- *)
   (* eta : ROUTE 1 -- MANUSCRIPT DEFINITION                    *)
   (* --------------------------------------------------------- *)

   etaFromDefinition =
      etaCFromDefinition*
      (dPhiDZValue^2)/
      (-GValue);


   (* --------------------------------------------------------- *)
   (* eta : ROUTE 2 -- EXPLICIT delta^4 CANCELLATION           *)
   (*                                                           *)
   (* eta = -C Fcc Lambda^2 G Omega^2                           *)
   (* --------------------------------------------------------- *)

   prefactorC =
      (sigmaInf mu^2 Jw l)/
      (4 PosDb epsilon^2 zeta^2);

   etaSimplified =
      -prefactorC*
       fccValue*
       lambdaValue^2*
       GValue*
       OmegaValue^2;


   etaAbsDifference =
      Abs[etaFromDefinition - etaSimplified];

   etaRelDifference =
      etaRelativeDifference[
         etaFromDefinition,
         etaSimplified
      ];


   (* --------------------------------------------------------- *)
   (* INTERNAL VERIFICATION                                     *)
   (* --------------------------------------------------------- *)

   verificationTolerance = 10^-12;

   etaCAgreementQ =
      TrueQ[etaCRelDifference <= verificationTolerance];

   etaAgreementQ =
      TrueQ[etaRelDifference <= verificationTolerance];

   If[
      ! etaAgreementQ,
      Message[
         ComputeEtaParametric::verify,
         N[etaRelDifference, 18]
      ]
   ];


   (* --------------------------------------------------------- *)
   (* RETURN OBJECT                                             *)
   (* --------------------------------------------------------- *)

   <|
      "Delta" -> deltaValue,
      "Lambda" -> lambdaValue,
      "PiD" -> piDValue,

      "Fcc" -> fccValue,
      "Omega" -> OmegaValue,
      "G" -> GValue,
      "dPhiDZ" -> dPhiDZValue,

      "phiC" -> phiCValue,

      "EtaC" -> etaCFromDefinition,
      "EtaCFromDefinition" -> etaCFromDefinition,
      "EtaCSimplified" -> etaCSimplified,
      "EtaCAbsoluteDifference" -> etaCAbsDifference,
      "EtaCRelativeDifference" -> etaCRelDifference,
      "EtaCAgreementQ" -> etaCAgreementQ,

      "Eta" -> etaFromDefinition,
      "EtaFromDefinition" -> etaFromDefinition,
      "EtaSimplified" -> etaSimplified,
      "EtaAbsoluteDifference" -> etaAbsDifference,
      "EtaRelativeDifference" -> etaRelDifference,
      "EtaAgreementQ" -> etaAgreementQ,

      "VerificationTolerance" -> verificationTolerance
   |>
 ];