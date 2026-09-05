(* ============================================================ *)
(* Code developed by Sánchez Lozano, G., Cruz-López C.-A., and  *)
(* F. Méndez, for the research:                                 *)
(* "Power generation in cylindrical microchannels with high     *)
(* surface zeta potential: a survey on viscoelectric effects    *)
(* over streaming potentials and efficiency"                    *)
(* Year: 2026                                                   *)
(* ------------------------------------------------------------ *)
(* GenerateDirectEtaParametricStudy.wl                          *)
(* General generator of eta(delta) studies using the DIRECT     *)
(* augmented solver                                             *)
(* ============================================================ *)

directEtaStudyDirectory = DirectoryName[$InputFileName];

If[
   directEtaStudyDirectory === "",
   Print[
      "ERROR: Run this file with Get[\"full_path\\GenerateDirectEtaParametricStudy.wl\"]."
   ];
   Abort[];
];


(* ============================================================ *)
(* LOAD THE VALIDATED DIRECT SOLVER                             *)
(* ============================================================ *)

directSolverCandidates = {
   FileNameJoin[{directEtaStudyDirectory, "DirectParametricSolver.wl"}],
   FileNameJoin[{directEtaStudyDirectory, "DirectParametricSolver(1).wl"}]
};

directSolverFile =
   SelectFirst[
      directSolverCandidates,
      FileExistsQ,
      Missing["NotFound"]
   ];

If[
   MissingQ[directSolverFile],
   Print[
      "ERROR: Could not find DirectParametricSolver.wl ",
      "(or DirectParametricSolver(1).wl) in ",
      directEtaStudyDirectory
   ];
   Abort[];
];

Get[directSolverFile];

Get[
   FileNameJoin[
      {directEtaStudyDirectory, "DirectEta_parametric.wl"}
   ]
];


ClearAll[
   GenerateDirectEtaParametricStudy,
   buildDirectEtaConfiguration,
   evaluateDirectEtaStudyPoint,
   directEtaParameterLabel,
   directEtaExactNumberQ
];


(* ============================================================ *)
(* MESSAGES                                                     *)
(* ============================================================ *)

GenerateDirectEtaParametricStudy::param =
   "Unknown parameter \"`1`\". Allowed parameters are \"PsiS\", \"omega\", \"Lambda\", and \"PiD\".";

GenerateDirectEtaParametricStudy::vals =
   "The sweep vector for \"`1`\" contains non-exact or non-numeric values: `2`.";

GenerateDirectEtaParametricStudy::delta =
   "deltaValues must contain positive exact numeric values. Invalid entries: `1`.";

GenerateDirectEtaParametricStudy::empty =
   "The parameter vector and deltaValues must both be non-empty.";

GenerateDirectEtaParametricStudy::verify =
   "Internal eta verification failed at delta = `1`, `2` = `3`.";

GenerateDirectEtaParametricStudy::eta =
   "Eta was not numeric at delta = `1`, `2` = `3`.";

GenerateDirectEtaParametricStudy::solver =
   "The direct solver failed to return an Association at delta = `1`, `2` = `3`.";

GenerateDirectEtaParametricStudy::etamodule =
   "ComputeDirectEtaParametric failed to return an Association at delta = `1`, `2` = `3`.";


(* ============================================================ *)
(* EXACT INPUT POLICY                                           *)
(* ============================================================ *)

directEtaExactNumberQ[x_] :=
   NumericQ[x] && Precision[x] === Infinity;


(* ============================================================ *)
(* CONFIGURATION BUILDER                                        *)
(* ============================================================ *)

buildDirectEtaConfiguration[
   parameterName_,
   parameterValue_
   ] :=
 Switch[
   parameterName,
   "PsiS",   <|"PsiS" -> parameterValue|>,
   "omega",  <|"omega" -> parameterValue|>,
   "Lambda", <|"Lambda" -> parameterValue|>,
   "PiD",    <|"PiD" -> parameterValue|>,
   _,        $Failed
 ];


(* ============================================================ *)
(* SINGLE TABLE POINT                                           *)
(* ============================================================ *)

evaluateDirectEtaStudyPoint[
   d_,
   parameterName_,
   parameterValue_,
   nGrid_ : 1000
   ] :=
 Module[
   {config, sol, etaSol, etaValue},

   (* Analytically regular points: do not call the direct solver. *)
   If[
      parameterName === "PiD" && TrueQ[parameterValue == 0],
      Return[0]
   ];

   If[
      parameterName === "PsiS" && TrueQ[parameterValue == 0],
      Return[0]
   ];


   config =
      buildDirectEtaConfiguration[
         parameterName,
         parameterValue
      ];

   If[
      config === $Failed,
      Return[
         Missing[
            "Failed",
            <|
               "Delta" -> d,
               "Parameter" -> parameterName,
               "ParameterValue" -> parameterValue,
               "Reason" -> "Configuration could not be built"
            |>
         ]
      ]
   ];


   sol =
      DirectParametricSolver[
         d,
         config,
         nGrid
      ];


   If[
      ! AssociationQ[sol],
      Message[
         GenerateDirectEtaParametricStudy::solver,
         d,
         parameterName,
         parameterValue
      ];
      Return[
         Missing[
            "Failed",
            <|
               "Delta" -> d,
               "Parameter" -> parameterName,
               "ParameterValue" -> parameterValue,
               "Reason" -> "Direct solver did not return an Association"
            |>
         ]
      ]
   ];


   etaSol =
      ComputeDirectEtaParametric[
         sol
      ];


   If[
      ! AssociationQ[etaSol],
      Message[
         GenerateDirectEtaParametricStudy::etamodule,
         d,
         parameterName,
         parameterValue
      ];
      Return[
         Missing[
            "Failed",
            <|
               "Delta" -> d,
               "Parameter" -> parameterName,
               "ParameterValue" -> parameterValue,
               "Reason" -> "ComputeDirectEtaParametric did not return an Association"
            |>
         ]
      ]
   ];


   If[
      ! (
         TrueQ[etaSol["EtaAgreementQ"]] &&
         TrueQ[etaSol["EtaCAgreementQ"]]
      ),
      Message[
         GenerateDirectEtaParametricStudy::verify,
         d,
         parameterName,
         parameterValue
      ];
      Return[
         Missing[
            "VerificationFailed",
            <|
               "Delta" -> d,
               "Parameter" -> parameterName,
               "ParameterValue" -> parameterValue,
               "EtaAgreementQ" -> etaSol["EtaAgreementQ"],
               "EtaCAgreementQ" -> etaSol["EtaCAgreementQ"],
               "EtaRelativeDifference" -> etaSol["EtaRelativeDifference"],
               "EtaCRelativeDifference" -> etaSol["EtaCRelativeDifference"]
            |>
         ]
      ]
   ];


   etaValue = etaSol["Eta"];


   If[
      NumericQ[etaValue],
      etaValue,
      Message[
         GenerateDirectEtaParametricStudy::eta,
         d,
         parameterName,
         parameterValue
      ];
      Missing[
         "Failed",
         <|
            "Delta" -> d,
            "Parameter" -> parameterName,
            "ParameterValue" -> parameterValue,
            "Reason" -> "Eta is not numeric"
         |>
      ]
   ]
 ];


(* ============================================================ *)
(* CLEAN COLUMN LABELS                                          *)
(* ============================================================ *)

directEtaParameterLabel[
   parameterName_,
   value_
   ] :=
   ToString[parameterName] <> "=" <> ToString[N[value]];


(* ============================================================ *)
(* GENERAL DIRECT ETA STUDY                                     *)
(* ============================================================ *)

GenerateDirectEtaParametricStudy[
   parameterName_,
   parameterValues_,
   deltaValues_,
   nGrid_ : 1000
   ] :=
 Module[
   {
      allowedParameters,
      invalidParameterValues,
      invalidDeltaValues,

      etaMatrix,
      failedPositions,
      labels,
      header,
      numericRows,
      exportTable,

      safeParameterName,
      outputFile,
      metadataFile,
      metadata,

      analyticPiDZeroPoints,
      analyticPsiSZeroPoints,
      successfulPoints,
      totalPoints,
      counter
   },


   allowedParameters =
      {"PsiS", "omega", "Lambda", "PiD"};


   If[
      ! MemberQ[allowedParameters, parameterName],
      Message[
         GenerateDirectEtaParametricStudy::param,
         parameterName
      ];
      Return[$Failed]
   ];


   If[
      parameterValues === {} || deltaValues === {},
      Message[
         GenerateDirectEtaParametricStudy::empty
      ];
      Return[$Failed]
   ];


   invalidParameterValues =
      Select[
         parameterValues,
         ! directEtaExactNumberQ[#] &
      ];

   If[
      invalidParameterValues =!= {},
      Message[
         GenerateDirectEtaParametricStudy::vals,
         parameterName,
         invalidParameterValues
      ];
      Return[$Failed]
   ];


   invalidDeltaValues =
      Select[
         deltaValues,
         ! (directEtaExactNumberQ[#] && # > 0) &
      ];

   If[
      invalidDeltaValues =!= {},
      Message[
         GenerateDirectEtaParametricStudy::delta,
         invalidDeltaValues
      ];
      Return[$Failed]
   ];


   totalPoints =
      Length[deltaValues]*
      Length[parameterValues];


   analyticPiDZeroPoints =
      If[
         parameterName === "PiD" &&
         MemberQ[parameterValues, 0],
         Length[deltaValues],
         0
      ];


   analyticPsiSZeroPoints =
      If[
         parameterName === "PsiS" &&
         MemberQ[parameterValues, 0],
         Length[deltaValues],
         0
      ];


   Print[""];
   Print["============================================================"];
   Print["GENERATING DIRECT eta PARAMETRIC STUDY"];
   Print["============================================================"];
   Print["Method: direct augmented sparse system"];
   Print["Parameter varied: ", parameterName];
   Print["delta values: ", deltaValues];
   Print["Parameter values: ", parameterValues];
   Print["nGrid: ", nGrid];
   Print["Total table points: ", totalPoints];

   If[
      analyticPiDZeroPoints > 0,
      Print[
         "Analytic PiD = 0 points: ",
         analyticPiDZeroPoints,
         " (Eta = 0; direct solver is not called)"
      ]
   ];

   If[
      analyticPsiSZeroPoints > 0,
      Print[
         "Analytic PsiS = 0 points: ",
         analyticPsiSZeroPoints,
         " (Eta = 0; Poisson-Boltzmann solve is bypassed)"
      ]
   ];

   Print[""];


   counter = 0;


   etaMatrix =
      Table[
         Table[

            counter++;

            PrintTemporary[
               "[",
               counter,
               "/",
               totalPoints,
               "]  delta = ",
               N[d, 8],
               ",  ",
               parameterName,
               " = ",
               N[parameterValue, 10]
            ];


            evaluateDirectEtaStudyPoint[
               d,
               parameterName,
               parameterValue,
               nGrid
            ],

            {parameterValue, parameterValues}
         ],

         {d, deltaValues}
      ];


   failedPositions =
      Position[
         etaMatrix,
         _Missing
      ];


   successfulPoints =
      totalPoints -
      Length[failedPositions];


   If[
      failedPositions === {},
      Print[""];
      Print["All direct eta table points are valid."],
      Print[""];
      Print[
         "WARNING: ",
         Length[failedPositions],
         " point(s) returned Missing."
      ];
      Print[
         "Failed matrix positions: ",
         failedPositions
      ]
   ];


   labels =
      directEtaParameterLabel[
         parameterName,
         #
      ] & /@ parameterValues;


   header =
      Join[
         {"delta"},
         labels
      ];


   numericRows =
      MapThread[
         Join[
            {#1},
            Map[
               If[
                  MissingQ[#],
                  #,
                  N[#, 17]
               ] &,
               #2
            ]
         ] &,
         {deltaValues, etaMatrix}
      ];


   exportTable =
      Prepend[
         numericRows,
         header
      ];


   safeParameterName =
      StringReplace[
         ToString[parameterName],
         {
            "_" -> "",
            " " -> ""
         }
      ];


   outputFile =
      FileNameJoin[
         {
            directEtaStudyDirectory,
            "Direct_Eta_vs_delta_varying_" <>
            safeParameterName <>
            ".txt"
         }
      ];


   Export[
      outputFile,
      exportTable,
      "TSV"
   ];


   metadata =
      {
         {"Observable", "Eta"},
         {"Method", "Direct augmented sparse system"},
         {"VariedParameter", parameterName},
         {"nGrid", nGrid},
         {"DeltaValuesExact", ToString[deltaValues, InputForm]},
         {"ParameterValuesExact", ToString[parameterValues, InputForm]},
         {"NumberOfTablePoints", totalPoints},
         {"SuccessfulPoints", successfulPoints},
         {"FailedPoints", Length[failedPositions]},
         {"AnalyticPiDZeroPoints", analyticPiDZeroPoints},
         {
            "PiDZeroTreatment",
            If[
               parameterName === "PiD" &&
               MemberQ[parameterValues, 0],
               "Eta=0 assigned analytically; direct solver bypassed",
               "Not applicable"
            ]
         },
         {"AnalyticPsiSZeroPoints", analyticPsiSZeroPoints},
         {
            "PsiSZeroTreatment",
            If[
               parameterName === "PsiS" &&
               MemberQ[parameterValues, 0],
               "Eta=0 assigned analytically; PB solve bypassed",
               "Not applicable"
            ]
         },
         {
            "LambdaTreatment",
            "When Lambda is explicitly swept, Fcc is still computed from the PB/electrochemical state and Lambda is imposed independently in DirectParametricSolver."
         },
         {
            "InternalEtaVerification",
            "EtaAgreementQ and EtaCAgreementQ required True for every non-analytic accepted point"
         },
         {
            "ExactInputPolicy",
            "All explicit sweep values and delta values must be exact integers or rationals"
         }
      };


   metadataFile =
      FileNameJoin[
         {
            directEtaStudyDirectory,
            "Direct_Eta_vs_delta_varying_" <>
            safeParameterName <>
            "_metadata.txt"
         }
      ];


   Export[
      metadataFile,
      metadata,
      "TSV"
   ];


   Print[""];
   Print["Data exported to:"];
   Print[outputFile];

   Print[""];
   Print["Metadata exported to:"];
   Print[metadataFile];

   Print[""];
   Print[
      "Accepted points: ",
      successfulPoints,
      "/",
      totalPoints
   ];

   Print[""];
   Print["Direct eta parametric study finalized."];


   <|
      "Observable" -> "Eta",
      "Method" -> "Direct augmented sparse system",
      "Parameter" -> parameterName,
      "ParameterValues" -> parameterValues,
      "DeltaValues" -> deltaValues,
      "nGrid" -> nGrid,

      "Matrix" -> etaMatrix,
      "Header" -> header,
      "Table" -> exportTable,

      "FailedPositions" -> failedPositions,
      "NumberOfTablePoints" -> totalPoints,
      "SuccessfulPoints" -> successfulPoints,

      "AnalyticPiDZeroPoints" -> analyticPiDZeroPoints,
      "AnalyticPsiSZeroPoints" -> analyticPsiSZeroPoints,

      "OutputFile" -> outputFile,
      "MetadataFile" -> metadataFile
   |>
 ];


(* ============================================================ *)
(* EXACT SWEEP VECTORS                                          *)
(* ============================================================ *)

DirectEtaPsiSValues = {
   -13/10,
   -123249/100000,
   -6/5,
   -11/10,
   -1,
   -9/10,
   -4/5,
   -7/10,
   -3/5,
   -1/2
};

DirectEtaOmegaValues =
   Range[0, 50, 5]/10000;

DirectEtaLambdaValues =
   Range[15, 105, 15]/100;

DirectEtaPiDValues =
   Range[0, 350, 50];

DirectEtaDeltaValues =
   Range[1, 30];