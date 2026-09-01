(* ============================================================ *)
(* Code developed by Sánchez Lozano, G., Cruz-López C.-A., and  *)
(* F. Méndez, for the research:                                 *)
(* "Power generation in cylindrical microchannels with high     *)
(* surface zeta potential: a survey on viscoelectric effects    *)
(* over streaming potentials and efficiency"                    *)
(* Year: 2026                                                   *)
(* ------------------------------------------------------------ *)
(* GenerateDirectPhiAvParametricStudy.wl                        *)
(* Unified front end for the DIRECT augmented solver            *)
(* ============================================================ *)

Print[""];
Print["============================================================"];
Print["GenerateDirectPhiAvParametricStudy.wl"];
Print["============================================================"];
Print[""];

ClearAll[
    GenerateDirectPhiAvParametricStudy,
    ExportDirectPhiAvStudy
];

studyDirectory = DirectoryName[$InputFileName];

Get[
    FileNameJoin[
        {
            studyDirectory,
            "DirectParametricSolver.wl"
        }
    ]
];


(* ------------------------------------------------------------ *)
(* Canonical grids                                              *)
(* ------------------------------------------------------------ *)

PhiAvDeltaValues =
    Range[1, 30];

PhiAvPsiSValues =
    {
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


(* ------------------------------------------------------------ *)
(* General direct study                                         *)
(* ------------------------------------------------------------ *)

GenerateDirectPhiAvParametricStudy[
    parameterName_,
    parameterValues_,
    deltaValues_ : PhiAvDeltaValues,
    nGrid_ : 1000
] := Module[
    {
        allowedParameters,
        rows,
        failures,
        result,
        config,
        phiAv,
        total,
        counter,
        safeParameterName,
        studyResult
    },

    allowedParameters =
        {"PsiS", "omega", "Lambda", "PiD"};

    If[
        !MemberQ[allowedParameters, parameterName],
        Print[
            "ERROR: parameterName must be one of ",
            allowedParameters,
            "."
        ];
        Return[$Failed];
    ];


    If[
        parameterName === "PiD",
        Print[
            "Study skipped: PhiAv is analytically independent of PiD."
        ];
        Return[
            <|
                "Status" -> "Skipped",
                "ParameterName" -> "PiD",
                "Reason" ->
                    "PhiAv is analytically independent of PiD."
            |>
        ];
    ];


    If[
        !And @@ (NumericQ /@ parameterValues) ||
        !And @@ (NumericQ /@ deltaValues),

        Print[
            "ERROR: parameterValues and deltaValues must be numeric."
        ];
        Return[$Failed];
    ];


    failures = {};
    total = Length[parameterValues]*Length[deltaValues];
    counter = 0;


    rows =
        Table[

            Prepend[

                Table[

                    counter++;

                    PrintTemporary[
                        "[",
                        counter,
                        "/",
                        total,
                        "]  delta = ",
                        N[delta, 8],
                        ",  ",
                        parameterName,
                        " = ",
                        N[param, 10]
                    ];


                    config =
                        <|
                            parameterName -> param
                        |>;


                    result =
                        DirectParametricSolver[
                            delta,
                            config,
                            nGrid
                        ];


                    phiAv =
                        If[
                            AssociationQ[result] &&
                            KeyExistsQ[result, "PhiAv"] &&
                            NumericQ[result["PhiAv"]],

                            result["PhiAv"],

                            AppendTo[
                                failures,
                                <|
                                    "Delta" -> delta,
                                    "ParameterName" -> parameterName,
                                    "ParameterValue" -> param
                                |>
                            ];

                            Missing["Failed"]
                        ];


                    phiAv,

                    {param, parameterValues}
                ],

                delta
            ],

            {delta, deltaValues}
        ];

    safeParameterName = StringReplace[parameterName, {"_" -> "", " " -> ""}];

    studyResult = <|
        "Status" -> "Completed",
        "Method" -> "Direct augmented sparse system",
        "ParameterName" -> parameterName,
        "ParameterValues" -> parameterValues,
        "DeltaValues" -> deltaValues,
        "nGrid" -> nGrid,
        "DataRows" -> rows,
        "Failures" -> failures
    |>;

    (* Automate export call within the generator *)
    ExportDirectPhiAvStudy[
        studyResult,
        "Direct_PhiAv_vs_delta_varying_" <> safeParameterName <> ".txt"
    ];

    studyResult
];


(* --------------------------------------------------------------------- *)
(* Export TSV in the same row/column structure used by analytic method   *)
(* --------------------------------------------------------------------- *)

ExportDirectPhiAvStudy[
    study_Association,
    fileName_String
] := Module[
    {
        parameterName,
        parameterValues,
        header,
        table,
        exportPath,
        metadataPath,
        metadata
    },

    If[
        Lookup[study, "Status", ""] =!= "Completed",
        Print["Nothing exported: the study was not completed."];
        Return[$Failed];
    ];


    parameterName = study["ParameterName"];
    parameterValues = study["ParameterValues"];

    header =
        Prepend[
            (
                parameterName <> "=" <> ToString[N[#]]
            ) & /@ parameterValues,
            "delta"
        ];

    table =
        Prepend[
            study["DataRows"],
            header
        ];


    exportPath =
        FileNameJoin[
            {studyDirectory, fileName}
        ];


    Export[
        exportPath,
        table,
        "TSV"
    ];


    metadataPath =
        FileNameJoin[
            {
                studyDirectory,
                FileBaseName[fileName] <> "_metadata.txt"
            }
        ];


    metadata =
        {
            "Method: Direct augmented sparse system",
            "Observable: Phi_av",
            "Parameter varied: " <> parameterName,
            "nGrid: " <> ToString[study["nGrid"]],
            "Delta values: " <> ToString[study["DeltaValues"], InputForm],
            "Parameter values: " <>
                ToString[study["ParameterValues"], InputForm],
            "Failures: " <>
                ToString[study["Failures"], InputForm]
        };


    Export[
        metadataPath,
        StringRiffle[metadata, "\n"],
        "Text"
    ];


    Print["Exported table: ", exportPath];
    Print["Exported metadata: ", metadataPath];

    If[
        Length[study["Failures"]] == 0,
        Print["No genuine solver failures were detected."],
        Print[
            "Failures detected: ",
            Length[study["Failures"]]
        ];
        Print[study["Failures"]];
    ];


    exportPath
];
