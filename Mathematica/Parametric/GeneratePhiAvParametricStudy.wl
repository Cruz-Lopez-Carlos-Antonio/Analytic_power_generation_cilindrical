(* ============================================================ *)
(* Code developed by Sánchez Lozano, G., Cruz-López C.-A., and  *)
(* F. Méndez, for the research:                                 *)
(* "Power generation in cylindrical microchannels with high     *)
(* surface zeta potential: a survey on viscoelectric effects    *)
(* over streaming potentials and efficiency"                    *)
(* Year: 2026                                                   *)
(* ------------------------------------------------------------ *)
(* GeneratePhiAvParametricStudy.wl                              *)
(* General generator of Phi_av(delta) parametric studies        *)
(* ============================================================ *)

studyDirectory = DirectoryName[$InputFileName];

If[
   studyDirectory === "",
   Print[
      "ERROR: Run this file with Get[\"full_path\\GeneratePhiAvParametricStudy.wl\"]."
   ];
   Abort[];
];

Get[
   FileNameJoin[
      {
         studyDirectory,
         "SemianalyticalParametricSolver_exact_inputs.wl"
      }
   ]
];

ClearAll[
   GeneratePhiAvParametricStudy,
   buildPhiAvConfiguration,
   evaluatePhiAvStudyPoint,
   parameterLabel
];

GeneratePhiAvParametricStudy::param =
   "Unknown parameter \"`1`\". Allowed parameters are \"PsiS\", \"omega\", and \"Lambda\".";

GeneratePhiAvParametricStudy::vals =
   "The sweep vector for \"`1`\" contains non-exact or non-numeric values: `2`.";

GeneratePhiAvParametricStudy::delta =
   "deltaValues must contain positive exact numeric values. Invalid entries: `1`.";

GeneratePhiAvParametricStudy::empty =
   "The parameter vector and deltaValues must both be non-empty.";

GeneratePhiAvParametricStudy::pid =
   "PiD sweep skipped: PhiAv is analytically independent of PiD in the present semianalytical formulation. The PiD dependence cancels exactly before evaluating PhiAv, so the corresponding curves are coincident. No numerical sweep will be launched.";

buildPhiAvConfiguration[parameterName_String, parameterValue_] :=
 Switch[
   parameterName,
   "PsiS",   <|"PsiS" -> parameterValue|>,
   "omega",  <|"omega" -> parameterValue|>,
   "Lambda", <|"Lambda" -> parameterValue|>,
   _,        $Failed
 ];

evaluatePhiAvStudyPoint[d_, parameterName_String, parameterValue_] :=
 Module[
   {config, result, value},

   config = buildPhiAvConfiguration[parameterName, parameterValue];

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

   result =
      SemianalyticalParametricSolver[
         d,
         config
      ];

   If[
      ! AssociationQ[result],
      Return[
         Missing[
            "Failed",
            <|
               "Delta" -> d,
               "Parameter" -> parameterName,
               "ParameterValue" -> parameterValue,
               "Reason" -> "Solver did not return an Association"
            |>
         ]
      ]
   ];

   If[
      ! KeyExistsQ[result, "PhiAv"],
      Return[
         Missing[
            "Failed",
            <|
               "Delta" -> d,
               "Parameter" -> parameterName,
               "ParameterValue" -> parameterValue,
               "Reason" -> "PhiAv key absent"
            |>
         ]
      ]
   ];

   value = result["PhiAv"];

   If[
      NumericQ[value],
      value,
      Missing[
         "Failed",
         <|
            "Delta" -> d,
            "Parameter" -> parameterName,
            "ParameterValue" -> parameterValue,
            "Reason" -> "PhiAv is not numeric"
         |>
      ]
   ]
 ];

(* Clean labels without NumberForm *)
parameterLabel[parameterName_String, value_] :=
   parameterName <> "=" <> ToString[N[value]];

GeneratePhiAvParametricStudy[
   parameterName_String,
   parameterValues_List,
   deltaValues_List
   ] :=
 Module[
   {
      allowedParameters,
      invalidParameterValues,
      invalidDeltaValues,
      phiAvMatrix,
      failedPositions,
      labels,
      header,
      numericRows,
      exportTable,
      safeParameterName,
      outputFile,
      metadataFile,
      metadata
   },

   allowedParameters = {"PsiS", "omega", "Lambda"};

   (* --------------------------------------------------------- *)
   (* PHYSICAL-LOGIC PROTECTION: PiD                            *)
   (* --------------------------------------------------------- *)

   If[
      parameterName === "PiD",
      Message[GeneratePhiAvParametricStudy::pid];
      Return[
         <|
            "Observable" -> "PhiAv",
            "Parameter" -> "PiD",
            "Status" -> "Skipped",
            "Reason" -> "PhiAv is analytically independent of PiD",
            "AnalyticalIndependence" -> True
         |>
      ]
   ];

   If[
      ! MemberQ[allowedParameters, parameterName],
      Message[GeneratePhiAvParametricStudy::param, parameterName];
      Return[$Failed]
   ];

   If[
      parameterValues === {} || deltaValues === {},
      Message[GeneratePhiAvParametricStudy::empty];
      Return[$Failed]
   ];

   invalidParameterValues =
      Select[
         parameterValues,
         !(NumericQ[#] && Precision[#] === Infinity) &
      ];

   If[
      invalidParameterValues =!= {},
      Message[
         GeneratePhiAvParametricStudy::vals,
         parameterName,
         invalidParameterValues
      ];
      Return[$Failed]
   ];

   invalidDeltaValues =
      Select[
         deltaValues,
         !(NumericQ[#] && Precision[#] === Infinity && # > 0) &
      ];

   If[
      invalidDeltaValues =!= {},
      Message[
         GeneratePhiAvParametricStudy::delta,
         invalidDeltaValues
      ];
      Return[$Failed]
   ];

   Print[""];
   Print["============================================================"];
   Print["GENERATING Phi_av PARAMETRIC STUDY"];
   Print["============================================================"];
   Print["Parameter varied: ", parameterName];
   Print["delta values: ", deltaValues];
   Print["Parameter values: ", parameterValues];
   Print[
      "Total evaluations: ",
      Length[deltaValues] Length[parameterValues]
   ];
   Print[""];

   phiAvMatrix =
      Table[
         Print[
            "delta = ",
            d,
            "   (",
            First@First@Position[deltaValues, d],
            "/",
            Length[deltaValues],
            ")"
         ];

         Table[
            evaluatePhiAvStudyPoint[
               d,
               parameterName,
               parameterValue
            ],
            {parameterValue, parameterValues}
         ],
         {d, deltaValues}
      ];

   failedPositions =
      Position[
         phiAvMatrix,
         _Missing
      ];

   If[
      failedPositions === {},
      Print[""];
      Print["All points returned numeric PhiAv values."],
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
      parameterLabel[
         parameterName,
         #
      ] & /@
      parameterValues;

   header =
      Join[
         {"delta"},
         labels
      ];

   numericRows =
      MapThread[
         Join[
            {#1},
            N[#2, 17]
         ] &,
         {
            deltaValues,
            phiAvMatrix
         }
      ];

   exportTable =
      Prepend[
         numericRows,
         header
      ];

   safeParameterName =
      StringReplace[
         parameterName,
         {
            "_" -> "",
            " " -> ""
         }
      ];

   outputFile =
      FileNameJoin[
         {
            studyDirectory,
            "PhiAv_vs_delta_varying_" <>
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
         {"Observable", "PhiAv"},
         {"VariedParameter", parameterName},
         {
            "DeltaValuesExact",
            ToString[deltaValues, InputForm]
         },
         {
            "ParameterValuesExact",
            ToString[parameterValues, InputForm]
         },
         {
            "NumberOfEvaluations",
            Length[deltaValues] Length[parameterValues]
         },
         {
            "FailedPoints",
            Length[failedPositions]
         }
      };

   metadataFile =
      FileNameJoin[
         {
            studyDirectory,
            "PhiAv_vs_delta_varying_" <>
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
   Print["Phi_av parametric study finalized."];

   <|
      "Observable" -> "PhiAv",
      "Parameter" -> parameterName,
      "ParameterValues" -> parameterValues,
      "DeltaValues" -> deltaValues,
      "Matrix" -> phiAvMatrix,
      "Header" -> header,
      "Table" -> exportTable,
      "FailedPositions" -> failedPositions,
      "OutputFile" -> outputFile,
      "MetadataFile" -> metadataFile
   |>
 ];

PhiAvPsiSValues = {
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

PhiAvOmegaValues =
   Range[0, 50, 5]/10000;

PhiAvLambdaValues =
   Range[15, 105, 15]/100;


PhiAvDeltaValues =
   Range[1, 30];
