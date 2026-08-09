(* ========================================================= *)

ClearAll[PsiPB, PsiPrimePB, pPB];


(* Dimensionless electrostatic potential *)

PsiPB[Rin_?NumericQ] :=
  Which[

    Rin == 0,
    centerPotential,

    0 < Rin <= 1,
    PsiFunction[Rin],

    True,
    Indeterminate
  ];


(* Auxiliary variable p(R) = R Psi'(R) *)

pPB[Rin_?NumericQ] :=
  Which[

    Rin == 0,
    0,

    0 < Rin <= 1,
    pFunction[Rin],

    True,
    Indeterminate
  ];


(* Radial derivative of the potential *)

PsiPrimePB[Rin_?NumericQ] :=
  Which[

    Rin == 0,
    0,

    0 < Rin <= 1,
    pFunction[Rin]/Rin,

    True,
    Indeterminate
  ];


