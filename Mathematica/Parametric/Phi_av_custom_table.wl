(* ============================================================ *)
(* Phi_av_custom_table_corrected.wl                              *)
(* Robust custom table for autonomous PhiAvParametric            *)
(* ============================================================ *)
(*
   IMPORTANT CORRECTION

   The previous version used:

       Quiet[Check[PhiAvParametric[d, ps], $Failed]]

   In Mathematica, Check returns its failure expression when ANY
   message is generated during the evaluation, even if the solver
   successfully returns a perfectly usable numerical result.

   Therefore benign numerical warnings from NIntegrate could be
   converted artificially into Missing["Failed"].

   This version:
      1. evaluates PhiAvParametric under Quiet;
      2. DOES NOT use Check to classify warnings as failures;
      3. accepts the result only if it is an Association containing
         a numerical "PhiAv" value;
      4. otherwise stores Missing["Failed"].
*)
(* ============================================================ *)

baseDir = DirectoryName[$InputFileName];

If[
   baseDir === "",
   Print[
      "ERROR: Run this file with Get[\"full_path\\Phi_av_custom_table_corrected.wl\"]."
   ];
   Abort[];
];

Get[
   FileNameJoin[
      {
         baseDir,
         "Phi_av_parametric.wl"
      }
   ]
];


(* ============================================================ *)
(* USER INPUT                                                    *)
(* ============================================================ *)

deltaValues = Range[30];

psiSValues = {
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


(* ============================================================ *)
(* INFORMATION                                                   *)
(* ============================================================ *)

Print[""];
Print["CUSTOM Phi_av PARAMETRIC TABLE"];
Print["delta values = ", deltaValues];
Print["PsiS values  = ", N[psiSValues, 8]];
Print[
   "Total calculations = ",
   Length[deltaValues] Length[psiSValues]
];
Print[""];


(* ============================================================ *)
(* COMPUTE MATRIX                                                *)
(* ============================================================ *)

phiAvMatrix =
   Table[

      Print[
         "delta = ",
         d,
         "   PsiS = ",
         N[ps, 8]
      ];


      result =
         Quiet[
            PhiAvParametric[
               d,
               ps
            ]
         ];


      If[
         AssociationQ[result] &&
         KeyExistsQ[result, "PhiAv"] &&
         NumericQ[result["PhiAv"]],

         result["PhiAv"],

         Missing["Failed"]
      ],

      {d, deltaValues},
      {ps, psiSValues}
   ];


phiAvMatrixNumeric =
   N[
      phiAvMatrix,
      12
   ];


(* ============================================================ *)
(* HEADERS                                                       *)
(* ============================================================ *)

headers =
   Prepend[
      (
         "PsiS=" <>
         ToString[
            NumberForm[
               N[#, 8],
               {10, 5}
            ],
            StandardForm
         ]
      ) & /@ psiSValues,
      "delta"
   ];


(* ============================================================ *)
(* COMPLETE NUMERICAL TABLE                                      *)
(* ============================================================ *)

tableRows =
   MapThread[
      Prepend,
      {
         phiAvMatrixNumeric,
         deltaValues
      }
   ];


completeTable =
   Prepend[
      tableRows,
      headers
   ];


(* ============================================================ *)
(* FORMATTED TABLE FOR NOTEBOOK                                  *)
(* ============================================================ *)

displayMatrix =
   Map[
      If[
         NumericQ[#],
         NumberForm[
            #,
            {14, 8}
         ],
         #
      ] &,
      phiAvMatrixNumeric,
      {2}
   ];


displayRows =
   MapThread[
      Prepend,
      {
         displayMatrix,
         deltaValues
      }
   ];


displayTable =
   Prepend[
      displayRows,
      headers
   ];


Print[""];

Print[
   Grid[
      displayTable,
      Frame -> All,
      Alignment -> Center,
      Spacings -> {1.2, 0.6}
   ]
];


(* ============================================================ *)
(* EXPORT                                                        *)
(* ============================================================ *)

outputFile =
   FileNameJoin[
      {
         baseDir,
         "Phi_av_table.txt"
      }
   ];


Export[
   outputFile,
   completeTable,
   "Table"
];


(* ============================================================ *)
(* SUMMARY                                                       *)
(* ============================================================ *)

failedPoints =
   Flatten[
      Table[
         If[
            MissingQ[phiAvMatrix[[i, j]]],
            {
               deltaValues[[i]],
               psiSValues[[j]]
            },
            Nothing
         ],
         {i, Length[deltaValues]},
         {j, Length[psiSValues]}
      ],
      1
   ];


Print[""];
Print["Calculations finished."];
Print["TXT table exported to:"];
Print[outputFile];

Print[""];

If[
   failedPoints === {},
   Print["No failed points were detected."],
   Print[
      "Failed points = ",
      failedPoints
   ]
];

Print[""];
Print["Raw matrix:       phiAvMatrix"];
Print["Numerical matrix: phiAvMatrixNumeric"];
Print["Complete table:   completeTable"];
