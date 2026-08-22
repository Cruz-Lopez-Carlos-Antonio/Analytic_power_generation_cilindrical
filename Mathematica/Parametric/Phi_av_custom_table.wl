(* ============================================================ *)
(* Phi_av_custom_table.txt                                      *)
(* ============================================================ *)
(* Place this file in the SAME folder as Phi_av_parametric.wl   *)
(* and its dependent parametric modules.                        *)
(* Run with Get["...\\Phi_av_custom_table.txt"].               *)
(* Edit only deltaValues and psiSValues below.                  *)
(* ============================================================ *)

baseDir = DirectoryName[$InputFileName];

If[baseDir === "",
   Print["ERROR: Run this file with Get[\"full_path\\Phi_av_custom_table.txt\"]."];
   Abort[];
];

Get[FileNameJoin[{baseDir, "Phi_av_parametric.wl"}]];

(* ================= USER INPUT ================= *)

deltaValues =Range[10]

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

(* ============================================== *)

Print[""];
Print["CUSTOM Phi_av PARAMETRIC TABLE"];
Print["delta values = ", deltaValues];
Print["PsiS values  = ", N[psiSValues, 8]];
Print["Total calculations = ", Length[deltaValues] Length[psiSValues]];
Print[""];

phiAvMatrix =
   Table[
      Print["delta = ", d, "   PsiS = ", N[ps, 8]];

      result = Quiet[
         Check[
            PhiAvParametric[d, ps],
            $Failed
         ]
      ];

      If[
         result === $Failed,
         Missing["Failed"],
         result["PhiAv"]
      ],

      {d, deltaValues},
      {ps, psiSValues}
   ];

phiAvMatrixNumeric = N[phiAvMatrix, 12];

headers =
   Prepend[
      ("PsiS=" <> ToString[NumberForm[N[#, 8], {10, 5}], InputForm]) & /@ psiSValues,
      "delta"
   ];

tableRows =
   MapThread[
      Prepend,
      {phiAvMatrixNumeric, deltaValues}
   ];

completeTable = Prepend[tableRows, headers];

displayMatrix =
   Map[
      If[NumericQ[#], NumberForm[#, {14, 8}], #] &,
      phiAvMatrixNumeric,
      {2}
   ];

displayRows =
   MapThread[
      Prepend,
      {displayMatrix, deltaValues}
   ];

displayTable = Prepend[displayRows, headers];

Print[""];
Print[
   Grid[
      displayTable,
      Frame -> All,
      Alignment -> Center,
      Spacings -> {1.2, 0.6}
   ]
];

outputFile = FileNameJoin[{baseDir, "Phi_av_table.txt"}];

Export[
   outputFile,
   completeTable,
   "Table"
];

Print[""];
Print["Calculations finished."];
Print["TXT table exported to:"];
Print[outputFile];
Print[""];
Print["Raw matrix:       phiAvMatrix"];
Print["Numerical matrix: phiAvMatrixNumeric"];
Print["Complete table:   completeTable"];
