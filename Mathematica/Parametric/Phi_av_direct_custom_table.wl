(* ============================================================ *)
(* Phi_av_direct_custom_table.wl                                *)
(* Robust table generator for Hipatia's direct solver           *)
(* ============================================================ *)
(*
   PURPOSE

   Builds a Phi_av(delta, PsiS) table using:

       DirectSolverParametric[delta, PsiS, nGrid]

   IMPORTANT CORRECTION

   This version does NOT use:

       Check[..., $Failed]

   around DirectSolverParametric.

   In Mathematica, Check can return its failure expression when
   the evaluated code emits a message, even if a usable numerical
   result is ultimately returned.

   Therefore, this version:
     1. evaluates the direct solver under Quiet;
     2. inspects the returned object;
     3. accepts it only if it is an Association containing a
        numerical value under the key "PhiAv";
     4. otherwise records Missing["Failed"].

   The code also prints the final table and exports it as:

       Phi_av_direct_table.txt

   in the same folder.
*)
(* ============================================================ *)


(* ============================================================ *)
(* 1. DIRECTORY                                                 *)
(* ============================================================ *)

baseDir = DirectoryName[$InputFileName];

If[
   baseDir === "",
   Print[
      "ERROR: execute this file using Get[\"full_path\\Phi_av_direct_custom_table.wl\"]."
   ];
   Abort[];
];


(* ============================================================ *)
(* 2. LOAD HIPATIA DIRECT-SOLVER DEPENDENCIES                   *)
(* ============================================================ *)

Get[
   FileNameJoin[
      {baseDir, "Parameters.wl"}
   ]
];

Get[
   FileNameJoin[
      {baseDir, "PoissonBoltzmann_parametric.wl"}
   ]
];

Get[
   FileNameJoin[
      {baseDir, "F_cc_parametric.wl"}
   ]
];

Get[
   FileNameJoin[
      {baseDir, "Lambda_parameter_parametric.wl"}
   ]
];

Get[
   FileNameJoin[
      {baseDir, "MR_parametric.wl"}
   ]
];

Get[
   FileNameJoin[
      {baseDir, "Direct_solver_parametric.wl"}
   ]
];


Print[""];
Print["Hipatia direct-solver modules loaded correctly."];


(* ============================================================ *)
(* 3. USER INPUT                                                *)
(* Edit only these three definitions.                           *)
(* ============================================================ *)

deltaValues = Range[10];

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

nGridTable = 1000;


(* ============================================================ *)
(* 4. INFORMATION                                               *)
(* ============================================================ *)

Print[""];
Print["============================================================"];
Print["         CUSTOM DIRECT-SOLVER Phi_av TABLE"];
Print["============================================================"];
Print[""];

Print["delta values = ", deltaValues];
Print["PsiS values  = ", N[psiSValues, 8]];
Print["nGrid        = ", nGridTable];

Print[""];
Print[
   "Total calculations = ",
   Length[deltaValues] Length[psiSValues]
];

Print[""];
Print["Starting direct-solver calculations..."];
Print[""];


(* ============================================================ *)
(* 5. SINGLE-POINT SAFE EVALUATOR                               *)
(* ============================================================ *)

ClearAll[DirectPhiAvValue];

DirectPhiAvValue[
   deltaValue_?NumericQ,
   psiSValue_?NumericQ
   ] :=
 Module[
   {
      sol,
      value
   },

   sol =
      Quiet[
         DirectSolverParametric[
            N[deltaValue, 15],
            N[psiSValue, 15],
            nGridTable
         ]
      ];

   If[
      AssociationQ[sol],

      value =
         Lookup[
            sol,
            "PhiAv",
            Missing["NoPhiAv"]
         ],

      Return[
         Missing["Failed"]
      ]
   ];

   If[
      NumericQ[value],
      value,
      Missing["Failed"]
   ]
];


(* ============================================================ *)
(* 6. COMPUTE Phi_av MATRIX                                     *)
(* ============================================================ *)

phiAvDirectMatrix =
   Table[

      Print[
         "delta = ",
         deltaValue,
         "   PsiS = ",
         N[psiSValue, 8]
      ];

      DirectPhiAvValue[
         deltaValue,
         psiSValue
      ],

      {deltaValue, deltaValues},
      {psiSValue, psiSValues}
   ];


(* ============================================================ *)
(* 7. NUMERICAL MATRIX                                          *)
(* ============================================================ *)

phiAvDirectMatrixNumeric =
   N[
      phiAvDirectMatrix,
      12
   ];


(* ============================================================ *)
(* 8. HEADERS                                                    *)
(* ============================================================ *)

headers =
   Prepend[
      Map[
         Function[
            ps,
            "PsiS=" <>
            ToString[
               NumberForm[
                  N[ps, 8],
                  {10, 5}
               ],
               StandardForm
            ]
         ],
         psiSValues
      ],
      "delta"
   ];


(* ============================================================ *)
(* 9. COMPLETE TABLE                                            *)
(* ============================================================ *)

tableRows =
   MapThread[
      Prepend,
      {
         phiAvDirectMatrixNumeric,
         deltaValues
      }
   ];

completeDirectTable =
   Prepend[
      tableRows,
      headers
   ];


(* ============================================================ *)
(* 10. FORMATTED NOTEBOOK TABLE                                 *)
(* ============================================================ *)

displayMatrix =
   Map[
      Function[
         x,
         If[
            NumericQ[x],
            NumberForm[x, {14, 8}],
            x
         ]
      ],
      phiAvDirectMatrixNumeric,
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

displayDirectTable =
   Prepend[
      displayRows,
      headers
   ];


Print[""];

Print[
   Grid[
      displayDirectTable,
      Frame -> All,
      Alignment -> Center,
      Spacings -> {1.2, 0.6}
   ]
];


(* ============================================================ *)
(* 11. IDENTIFY GENUINE FAILED POINTS                           *)
(* ============================================================ *)

failedDirectPoints =
   Flatten[
      Table[
         If[
            MissingQ[
               phiAvDirectMatrix[[i, j]]
            ],

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


(* ============================================================ *)
(* 12. EXPORT TXT TABLE                                         *)
(* ============================================================ *)

outputFile =
   FileNameJoin[
      {
         baseDir,
         "Phi_av_direct_table.txt"
      }
   ];

Export[
   outputFile,
   completeDirectTable,
   "Table"
];


(* ============================================================ *)
(* 13. FINAL SUMMARY                                            *)
(* ============================================================ *)

Print[""];
Print["============================================================"];
Print["Direct-solver calculations finished."];
Print[""];

Print["TXT table exported to:"];
Print[outputFile];

Print[""];

If[
   failedDirectPoints === {},

   Print[
      "No genuine failed points were detected."
   ],

   Print[
      "Genuine failed points = ",
      failedDirectPoints
   ]
];

Print[""];

Print["Raw matrix:       phiAvDirectMatrix"];
Print["Numerical matrix: phiAvDirectMatrixNumeric"];
Print["Complete table:   completeDirectTable"];
Print["Formatted table:  displayDirectTable"];

Print["============================================================"];
