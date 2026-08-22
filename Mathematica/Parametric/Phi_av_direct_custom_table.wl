(* ============================================================ *)
(* Phi_av_direct_custom_table.txt                               *)
(* Custom Phi_av table using Hipatia's direct parametric solver *)
(* ============================================================ *)
(*
   PURPOSE

   This file builds a Phi_av(delta, PsiS) table using the
   direct parametric route implemented through

       Direct_solver_parametric.wl

   It is intentionally separate from the semianalytical
   Phi_av_custom_table.txt workflow.

   HOW TO USE

   1. Place this file in the SAME folder as:

          Parameters.wl
          PoissonBoltzmann_parametric.wl
          F_cc_parametric.wl
          Lambda_parameter_parametric.wl
          MR_parametric.wl
          Direct_solver_parametric.wl

   2. Edit only:

          deltaValues
          psiSValues
          nGridTable

   3. Run, for example:

          Get["C:\\Users\\AMD RYZEN 7\\Downloads\\Phi_av_direct_custom_table.txt"]

   4. The routine:
        - evaluates every requested pair (delta, PsiS),
        - prints a formatted table in Mathematica,
        - exports a TXT table in the same folder.

   The structure follows the dependencies used in Hipatia's
   Table_1_final.txt.
*)
(* ============================================================ *)


(* ============================================================ *)
(* 1. DIRECTORY                                                 *)
(* ============================================================ *)

baseDir = DirectoryName[$InputFileName];

If[
   baseDir === "",
   Print[
      "ERROR: Run this file with Get[\"full_path\\Phi_av_direct_custom_table.txt\"]."
   ];
   Abort[];
];


(* ============================================================ *)
(* 2. LOAD HIPATIA DIRECT-SOLVER DEPENDENCIES                   *)
(* ============================================================ *)

Get[
   FileNameJoin[
      {
         baseDir,
         "Parameters.wl"
      }
   ]
];

Get[
   FileNameJoin[
      {
         baseDir,
         "PoissonBoltzmann_parametric.wl"
      }
   ]
];

Get[
   FileNameJoin[
      {
         baseDir,
         "F_cc_parametric.wl"
      }
   ]
];

Get[
   FileNameJoin[
      {
         baseDir,
         "Lambda_parameter_parametric.wl"
      }
   ]
];

Get[
   FileNameJoin[
      {
         baseDir,
         "MR_parametric.wl"
      }
   ]
];

Get[
   FileNameJoin[
      {
         baseDir,
         "Direct_solver_parametric.wl"
      }
   ]
];


Print[""];
Print["Hipatia direct-solver modules loaded correctly."];


(* ============================================================ *)
(* 3. USER INPUT                                                *)
(* ============================================================ *)
(* Edit these values only.                                      *)

deltaValues =  Range[10];


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


(* Grid size used by DirectSolverParametric.                     *)
(* Hipatia's complete-table routine used nGrid = 1000.          *)

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
(* 5. COMPUTE Phi_av MATRIX                                     *)
(* ============================================================ *)

phiAvDirectMatrix =
   Table[

      Print[
         "delta = ",
         deltaValue,
         "   PsiS = ",
         N[psiSValue, 8]
      ];


      currentDirectSolution =
         Quiet[
            Check[
               DirectSolverParametric[
                  N[deltaValue, 15],
                  N[psiSValue, 15],
                  nGridTable
               ],
               $Failed
            ]
         ];


      If[
         currentDirectSolution === $Failed,

         Missing["Failed"],

         currentDirectSolution["PhiAv"]
      ],

      {deltaValue, deltaValues},

      {psiSValue, psiSValues}
   ];


(* ============================================================ *)
(* 6. NUMERICAL MATRIX                                          *)
(* ============================================================ *)

phiAvDirectMatrixNumeric =
   N[
      phiAvDirectMatrix,
      12
   ];


(* ============================================================ *)
(* 7. HEADERS                                                    *)
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
            InputForm
         ]
      ) & /@ psiSValues,
      "delta"
   ];


(* ============================================================ *)
(* 8. BUILD COMPLETE TABLE                                      *)
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
(* 9. FORMATTED TABLE FOR MATHEMATICA                           *)
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
(* 10. EXPORT TXT TABLE                                         *)
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
(* 11. FINAL MESSAGE                                            *)
(* ============================================================ *)

Print[""];
Print["============================================================"];
Print["Direct-solver calculations finished."];
Print[""];
Print["TXT table exported to:"];
Print[outputFile];
Print[""];
Print["Raw matrix:       phiAvDirectMatrix"];
Print["Numerical matrix: phiAvDirectMatrixNumeric"];
Print["Complete table:   completeDirectTable"];
Print["Formatted table:  displayDirectTable"];
Print["============================================================"];
