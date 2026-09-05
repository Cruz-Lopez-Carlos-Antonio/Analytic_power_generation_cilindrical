(* ============================================================ *)
(* Code developed by Sánchez Lozano, G., Cruz-López C.-A., and  *)
(* F. Méndez, for the research:                                 *)
(* "Power generation in cylindrical microchannels with high     *)
(* surface zeta potential: a survey on viscoelectric effects    *)
(* over streaming potentials and efficiency"                    *)
(* Year: 2026                                                   *)
(* ------------------------------------------------------------ *)
(* GenerateEtaParametricStudy.wl                                *)
(* General generator of eta(delta) parametric studies           *)
(* ============================================================ *)

etaStudyDirectory = DirectoryName[$InputFileName];

If[
   etaStudyDirectory === "",
   Print["ERROR: Run this file with Get[\"full_path\\GenerateEtaParametricStudy.wl\"]."];
   Abort[];
];

Get[FileNameJoin[{etaStudyDirectory, "SemianalyticalParametricSolver_exact_inputs.wl"}]];
Get[FileNameJoin[{etaStudyDirectory, "Eta_parametric.wl"}]];

ClearAll[
   GenerateEtaParametricStudy,
   buildEtaConfiguration,
   evaluateEtaStudyPoint,
   etaParameterLabel,
   etaExactNumberQ
];

GenerateEtaParametricStudy::param = "Unknown parameter \"`1`\". Allowed parameters are \"PsiS\", \"omega\", \"Lambda\", and \"PiD\".";
GenerateEtaParametricStudy::vals = "The sweep vector for \"`1`\" contains non-exact or non-numeric values: `2`.";
GenerateEtaParametricStudy::delta = "deltaValues must contain positive exact numeric values. Invalid entries: `1`.";
GenerateEtaParametricStudy::empty = "The parameter vector and deltaValues must both be non-empty.";
GenerateEtaParametricStudy::verify = "Internal eta verification failed at delta = `1`, `2` = `3`.";
GenerateEtaParametricStudy::eta = "Eta was not numeric at delta = `1`, `2` = `3`.";
GenerateEtaParametricStudy::solver = "The semianalytical solver failed to return an Association at delta = `1`, `2` = `3`.";
GenerateEtaParametricStudy::etamodule = "ComputeEtaParametric failed to return an Association at delta = `1`, `2` = `3`.";

etaExactNumberQ[x_] := NumericQ[x] && Precision[x] === Infinity;

buildEtaConfiguration[parameterName_String, parameterValue_] :=
 Switch[
   parameterName,
   "PsiS",   <|"PsiS" -> parameterValue|>,
   "omega",  <|"omega" -> parameterValue|>,
   "Lambda", <|"Lambda" -> parameterValue|>,
   "PiD",    <|"PiD" -> parameterValue|>,
   _,        $Failed
 ];

evaluateEtaStudyPoint[d_, parameterName_String, parameterValue_] :=
 Module[
   {config, sol, etaSol, etaValue},

   If[parameterName === "PiD" && TrueQ[parameterValue == 0], Return[0]];
   If[parameterName === "PsiS" && TrueQ[parameterValue == 0], Return[0]];

   config = buildEtaConfiguration[parameterName, parameterValue];
   If[config === $Failed, Return[Missing["Failed", <|"Delta" -> d, "Parameter" -> parameterName, "ParameterValue" -> parameterValue, "Reason" -> "Configuration could not be built"|>]]];

   sol = SemianalyticalParametricSolver[d, config];
   If[! AssociationQ[sol],
      Message[GenerateEtaParametricStudy::solver, d, parameterName, parameterValue];
      Return[Missing["Failed", <|"Delta" -> d, "Parameter" -> parameterName, "ParameterValue" -> parameterValue, "Reason" -> "Solver did not return an Association"|>]]
   ];

   etaSol = ComputeEtaParametric[sol];
   If[! AssociationQ[etaSol],
      Message[GenerateEtaParametricStudy::etamodule, d, parameterName, parameterValue];
      Return[Missing["Failed", <|"Delta" -> d, "Parameter" -> parameterName, "ParameterValue" -> parameterValue, "Reason" -> "ComputeEtaParametric did not return an Association"|>]]
   ];

   If[! (TrueQ[etaSol["EtaAgreementQ"]] && TrueQ[etaSol["EtaCAgreementQ"]]),
      Message[GenerateEtaParametricStudy::verify, d, parameterName, parameterValue];
      Return[Missing["VerificationFailed", <|"Delta" -> d, "Parameter" -> parameterName, "ParameterValue" -> parameterValue, "EtaAgreementQ" -> etaSol["EtaAgreementQ"], "EtaCAgreementQ" -> etaSol["EtaCAgreementQ"], "EtaRelativeDifference" -> etaSol["EtaRelativeDifference"], "EtaCRelativeDifference" -> etaSol["EtaCRelativeDifference"]|>]]
   ];

   etaValue = etaSol["Eta"];
   If[NumericQ[etaValue],
      etaValue,
      Message[GenerateEtaParametricStudy::eta, d, parameterName, parameterValue];
      Missing["Failed", <|"Delta" -> d, "Parameter" -> parameterName, "ParameterValue" -> parameterValue, "Reason" -> "Eta is not numeric"|>]
   ]
 ];

(* Clean labels without NumberForm *)
etaParameterLabel[parameterName_String, value_] := parameterName <> "=" <> ToString[N[value]];

GenerateEtaParametricStudy[
   parameterName_,
   parameterValues_,
   deltaValues_
   ] :=
 Module[
   {
      allowedParameters, invalidParameterValues, invalidDeltaValues,
      etaMatrix, failedPositions, labels, header, numericRows, exportTable,
      safeParameterName, outputFile, metadataFile, metadata,
      analyticPiDZeroPoints, analyticPsiSZeroPoints, successfulPoints, totalPoints,
      counter
   },

   allowedParameters = {"PsiS", "omega", "Lambda", "PiD"};

   If[! MemberQ[allowedParameters, parameterName], Message[GenerateEtaParametricStudy::param, parameterName]; Return[$Failed]];
   If[parameterValues === {} || deltaValues === {}, Message[GenerateEtaParametricStudy::empty]; Return[$Failed]];

   invalidParameterValues = Select[parameterValues, ! etaExactNumberQ[#] &];
   If[invalidParameterValues =!= {}, Message[GenerateEtaParametricStudy::vals, parameterName, invalidParameterValues]; Return[$Failed]];

   invalidDeltaValues = Select[deltaValues, ! (etaExactNumberQ[#] && # > 0) &];
   If[invalidDeltaValues =!= {}, Message[GenerateEtaParametricStudy::delta, invalidDeltaValues]; Return[$Failed]];

   totalPoints = Length[deltaValues] * Length[parameterValues];
   analyticPiDZeroPoints = If[parameterName === "PiD" && MemberQ[parameterValues, 0], Length[deltaValues], 0];
   analyticPsiSZeroPoints = If[parameterName === "PsiS" && MemberQ[parameterValues, 0], Length[deltaValues], 0];

   Print[""];
   Print["============================================================"];
   Print["GENERATING eta PARAMETRIC STUDY"];
   Print["============================================================"];
   Print["Parameter varied: ", parameterName];
   Print["delta values: ", deltaValues];
   Print["Parameter values: ", parameterValues];
   Print["Total table points: ", totalPoints];
   If[analyticPiDZeroPoints > 0, Print["Analytic PiD = 0 points: ", analyticPiDZeroPoints, " (Eta = 0; general solver is not called)"]];
   If[analyticPsiSZeroPoints > 0, Print["Analytic PsiS = 0 points: ", analyticPsiSZeroPoints, " (Eta = 0; Poisson-Boltzmann FindRoot is bypassed)"]];
   Print[""];

   counter = 0;
   etaMatrix =
      Table[
         Table[
            counter++;
            PrintTemporary[
               "[", counter, "/", totalPoints, "]  delta = ",
               N[d, 8], ",  ", parameterName, " = ", N[parameterValue, 10]
            ];
            evaluateEtaStudyPoint[d, parameterName, parameterValue],
            {parameterValue, parameterValues}
         ],
         {d, deltaValues}
      ];

   failedPositions = Position[etaMatrix, _Missing];
   successfulPoints = totalPoints - Length[failedPositions];

   If[failedPositions === {},
      Print[""]; Print["All eta table points are valid."],
      Print[""]; Print["WARNING: ", Length[failedPositions], " point(s) returned Missing."];
      Print["Failed matrix positions: ", failedPositions]
   ];

   labels = etaParameterLabel[parameterName, #] & /@ parameterValues;
   header = Join[{"delta"}, labels];
   numericRows = MapThread[Join[{#1}, Map[If[MissingQ[#], #, N[#, 17]] &, #2]] &, {deltaValues, etaMatrix}];
   exportTable = Prepend[numericRows, header];

   safeParameterName = StringReplace[parameterName, {"_" -> "", " " -> ""}];
   outputFile = FileNameJoin[{etaStudyDirectory, "Eta_vs_delta_varying_" <> safeParameterName <> ".txt"}];
   Export[outputFile, exportTable, "TSV"];

   metadata = {
         {"Observable", "Eta"},
         {"VariedParameter", parameterName},
         {"DeltaValuesExact", ToString[deltaValues, InputForm]},
         {"ParameterValuesExact", ToString[parameterValues, InputForm]},
         {"NumberOfTablePoints", totalPoints},
         {"SuccessfulPoints", successfulPoints},
         {"FailedPoints", Length[failedPositions]},
         {"AnalyticPiDZeroPoints", analyticPiDZeroPoints},
         {"PiDZeroTreatment", If[parameterName === "PiD" && MemberQ[parameterValues, 0], "Eta=0 assigned analytically; general solver bypassed", "Not applicable"]},
         {"AnalyticPsiSZeroPoints", analyticPsiSZeroPoints},
         {"PsiSZeroTreatment", If[parameterName === "PsiS" && MemberQ[parameterValues, 0], "Eta=0 assigned analytically; Poisson-Boltzmann FindRoot bypassed", "Not applicable"]},
         {"InternalEtaVerification", "EtaAgreementQ and EtaCAgreementQ required True for every non-analytic accepted point"},
         {"ExactInputPolicy", "All explicit sweep values and delta values must be exact integers or rationals"}
      };

   metadataFile = FileNameJoin[{etaStudyDirectory, "Eta_vs_delta_varying_" <> safeParameterName <> "_metadata.txt"}];
   Export[metadataFile, metadata, "TSV"];

   Print[""]; Print["Data exported to:"]; Print[outputFile];
   Print[""]; Print["Metadata exported to:"]; Print[metadataFile];
   Print[""]; Print["Accepted points: ", successfulPoints, "/", totalPoints];
   Print[""]; Print["eta parametric study finalized."];

   <|
      "Observable" -> "Eta", "Parameter" -> parameterName, "ParameterValues" -> parameterValues, "DeltaValues" -> deltaValues,
      "Matrix" -> etaMatrix, "Header" -> header, "Table" -> exportTable, "FailedPositions" -> failedPositions,
      "NumberOfTablePoints" -> totalPoints, "SuccessfulPoints" -> successfulPoints,
      "AnalyticPiDZeroPoints" -> analyticPiDZeroPoints, "AnalyticPsiSZeroPoints" -> analyticPsiSZeroPoints,
      "OutputFile" -> outputFile, "MetadataFile" -> metadataFile
   |>
 ];

(* ============================================================ *)
(* EXACT SWEEP VECTORS   *)
(* ============================================================ *)

EtaPsiSValues = {-13/10, -123249/100000, -6/5, -11/10, -1, -9/10, -4/5, -7/10, -3/5, -1/2};

EtaOmegaValues = Range[0, 50, 5]/10000;

EtaLambdaValues = Range[15, 105, 15]/100;

EtaPiDValues = Range[50, 350, 50];

EtaDeltaValues = Range[1, 30];