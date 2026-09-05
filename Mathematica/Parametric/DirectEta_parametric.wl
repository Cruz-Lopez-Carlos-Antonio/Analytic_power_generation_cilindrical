(* ============================================================ *)
(* Code developed by Sánchez Lozano, G., Cruz-López C.-A., and  *)
(* F. Méndez, for the research:                                 *)
(* "Power generation in cylindrical microchannels with high     *)
(* surface zeta potential: a survey on viscoelectric effects    *)
(* over streaming potentials and efficiency"                    *)
(* Year: 2026                                                   *)
(* ------------------------------------------------------------ *)
(* DirectEta_parametric.wl                                      *)
(* Postprocessing of eta from the DIRECT augmented solver       *)
(* ============================================================ *)

ClearAll[
   ComputeDirectEtaParametric,
   directEtaRelativeDifference
];


(* ============================================================ *)
(* MESSAGES                                                     *)
(* ============================================================ *)

ComputeDirectEtaParametric::keys =
   "The supplied direct-solver Association is missing required key(s): `1`.";

ComputeDirectEtaParametric::phys =
   "One or more physical parameters required for eta are not numeric. Check Parameters.wl (sigmaInf, mu, Jw, l, epsilon, zeta, PosDb).";

ComputeDirectEtaParametric::pid0 =
   "PiD = 0 is a removable singular case. In a parametric study assign Eta -> 0 analytically before calling the direct solver.";

ComputeDirectEtaParametric::verify =
   "The two eta evaluations differ by more than the internal verification tolerance. Relative difference = `1`.";


(* ============================================================ *)
(* RELATIVE DIFFERENCE HELPER                                   *)
(* ============================================================ *)

directEtaRelativeDifference[x_, y_] :=
 Module[{scale},
   scale = Max[Abs[x], Abs[y]];
   If[TrueQ[scale == 0], 0, Abs[x - y]/scale]
 ];


(* ============================================================ *)
(* MAIN POSTPROCESSOR                                           *)
(* ============================================================ *)

ComputeDirectEtaParametric[sol_] :=
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

   If[! AssociationQ[sol],
      Return[$Failed]
   ];

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
         ComputeDirectEtaParametric::keys,
         missingKeys
      ];
      Return[$Failed]
   ];


   deltaValue  = sol["Delta"];
   lambdaValue = sol["Lambda"];
   piDValue    = sol["PiD"];
   fccValue    = sol["Fcc"];
   OmegaValue  = sol["Omega"];
   GValue      = sol["G"];
   dPhiDZValue = sol["dPhiDZ"];


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
      Message[ComputeDirectEtaParametric::phys];
      Return[$Failed]
   ];


   If[
      TrueQ[piDValue == 0],
      Message[ComputeDirectEtaParametric::pid0]
   ];


   (* Characteristic potential scale *)
   phiCValue =
      (mu Jw l)/
      (epsilon zeta deltaValue^2);


   (* eta_c, route 1: original definition *)
   etaCFromDefinition =
      ((1/4) sigmaInf fccValue (phiCValue^2/l^2))/
      (Jw (PosDb/l));


   (* eta_c, route 2: simplified form *)
   etaCSimplified =
      (sigmaInf fccValue mu^2 Jw l)/
      (4 PosDb epsilon^2 zeta^2 deltaValue^4);


   etaCAbsDifference =
      Abs[etaCFromDefinition - etaCSimplified];

   etaCRelDifference =
      directEtaRelativeDifference[
         etaCFromDefinition,
         etaCSimplified
      ];


   (* eta, route 1: manuscript definition *)
   etaFromDefinition =
      etaCFromDefinition*
      (dPhiDZValue^2)/
      (-GValue);


   (* eta, route 2: explicit delta^4 cancellation *)
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
      directEtaRelativeDifference[
         etaFromDefinition,
         etaSimplified
      ];


   verificationTolerance = 10^-12;

   etaCAgreementQ =
      TrueQ[etaCRelDifference <= verificationTolerance];

   etaAgreementQ =
      TrueQ[etaRelDifference <= verificationTolerance];


   If[
      ! etaAgreementQ,
      Message[
         ComputeDirectEtaParametric::verify,
         N[etaRelDifference, 18]
      ]
   ];


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