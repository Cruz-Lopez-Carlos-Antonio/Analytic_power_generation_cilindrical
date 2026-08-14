(* Corrected operational F.wl *)

fBaseDir = DirectoryName[$InputFileName];

Get[FileNameJoin[{fBaseDir, "F0.wl"}]];
Get[FileNameJoin[{fBaseDir, "F1.wl"}]];
Get[FileNameJoin[{fBaseDir, "Omega_parameter.wl"}]];

ClearAll[F, FPrime];

F[R_?NumericQ] := Module[{value},
  If[R < 0 || R > 1,
    Message[F::domain, R];
    Return[$Failed];
  ];

  value = F0[R] + Omega F1[R];

  If[R == 1, 0, value]
];

F::domain =
  "The argument R = `1` is outside the required domain 0 <= R <= 1.";

FPrime[R_?NumericQ] := Module[{value},
  If[R < 0 || R > 1,
    Message[FPrime::domain, R];
    Return[$Failed];
  ];

  value = F0Prime[R] + Omega F1Prime[R];

  If[R == 0, 0, value]
];

FPrime::domain =
  "The argument R = `1` is outside the required domain 0 <= R <= 1.";
