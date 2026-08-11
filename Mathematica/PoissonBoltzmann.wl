(* ========================================================= *)
(* POISSON-BOLTZMANN SOLVER                                  *)
(* ========================================================= *)
(*
   Cylindrical nonlinear Poisson-Boltzmann equation:
       Psi''(R) + (1/R) Psi'(R)= delta^2 Sinh[Psi(R)]
   0 <= R <= 1
   
   Boundary conditions:
       Psi'(0) = 0
       Psi(1)  = PsiS

   We introduce
       p(R) = R Psi'(R)

   so that
       Psi'(R) = p(R)/R
       p'(R)   = delta^2 R Sinh[Psi(R)]

   The apparent singularity at R = 0 is handled by
   starting at a small positive epsPB and using the
   regular Taylor expansion about the symmetry axis.
*)
(* ========================================================= *)


(* ========================================================= *)
(* 1. LOAD PARAMETERS                                        *)
(* ========================================================= *)

moduleDirectory = DirectoryName[$InputFileName];

Get[
  FileNameJoin[
    {moduleDirectory, "Parameters.wl"}
  ]
];

(* ========================================================= *)
(* 2. NUMERICAL PRECISION                                    *)
(* ========================================================= *)

(*
   The shooting problem is sensitive because Psi(0) is of
   order 10^-9 while Psi(1) is of order unity.

   Therefore, internal calculations are performed with
   higher working precision.
*)

wpPB = 30;
wpExt = 40;  (* Extended precision used to absorb numerical precision loss *)

deltaPB = SetPrecision[delta, wpExt];
PsiSPB  = SetPrecision[PsiS, wpExt];

(* Small distance from the symmetry axis *)
epsPB = SetPrecision[10^-8, wpExt];

(* ========================================================= *)
(* 3. LOCAL REGULAR EXPANSION AT R = 0                       *)
(* ========================================================= *)

(*
   Let Psi(0) = psi0.
   Regularity gives
       Psi(R)= psi0+ delta^2 Sinh[psi0] R^2/4+ O(R^4),

   and therefore

       p(R) = R Psi'(R)
            = delta^2 Sinh[psi0] R^2/2
              + O(R^4).

   These expressions supply the initial conditions at R = epsPB.
*)

(* ========================================================= *)
(* 4. INTEGRATION FOR A GIVEN CENTER POTENTIAL               *)
(* ========================================================= *)

ClearAll[pbForCenter, psi, p, R];
pbForCenter[psi0_?NumericQ] :=
  Module[
    {psi0PB},
    (* Force every shooting value to the solver precision using extended precision *)
    psi0PB = SetPrecision[psi0, wpExt];
    NDSolveValue[
      {
        psi'[R] == p[R]/R,
        p'[R] ==
          deltaPB^2*R*Sinh[psi[R]],
          
        (* Regular initial value of Psi near the axis *)
        psi[epsPB] ==psi0PB
          +
          (deltaPB^2*Sinh[psi0PB]/4)*epsPB^2,

        (* Regular initial value of p = R Psi' *)

        p[epsPB] == (deltaPB^2*Sinh[psi0PB]/2)*epsPB^2

      },
      {psi, p},
      {R, epsPB, 1},
      WorkingPrecision -> wpPB,
      AccuracyGoal -> 18,
      PrecisionGoal -> 18,
      MaxSteps -> Infinity
    ]
  ];

(* ========================================================= *)
(* 5. SHOOTING RESIDUAL                                      *)
(* ========================================================= *)

(*For an assumed center potential psi0, integrate
   to R = 1 and measure the mismatch residual(psi0)
       = Psi(1; psi0) - PsiS.
   The correct center potential makes this quantity zero.
*)

ClearAll[wallResidual];
wallResidual[psi0_?NumericQ] :=
  Module[{sol},sol = pbForCenter[psi0];
    sol[[1]][1] - PsiSPB];

(* ========================================================= *)
(* 6. INITIAL ESTIMATE FROM LINEARIZED DEBYE-HUCKEL THEORY   *)
(* ========================================================= *)

(*
   The linearized cylindrical solution is
       Psi(R) =PsiS I0(delta R)/I0(delta).
   Hence
       Psi(0) = PsiS/I0(delta).
   This is used ONLY as an initial numerical estimate.
*)

psi0Guess =
  PsiSPB/BesselI[0, deltaPB];

(* ========================================================= *)
(* 7. DETERMINE THE CENTER POTENTIAL                         *)
(* ========================================================= *)

ClearAll[psi0];
centerPotential =
  psi0 /. FindRoot[
    wallResidual[psi0] == 0,
    {
      psi0,
      (9/10)*psi0Guess,
      (11/10)*psi0Guess
    },
    WorkingPrecision -> wpPB,
    AccuracyGoal -> 20,
    PrecisionGoal -> 20,
    MaxIterations -> 100
  ];


(* ========================================================= *)
(* 8. FINAL POISSON-BOLTZMANN SOLUTION                       *)
(* ========================================================= *)

pbSolution =
  pbForCenter[centerPotential];
PsiFunction = pbSolution[[1]];
pFunction = pbSolution[[2]];

(* ========================================================= *)
(* 9. FUNCTIONS INCLUDING THE AXIS R = 0                     *)
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
    Rin == 0,0,
    0 < Rin <= 1,
    pFunction[Rin]/Rin,
    True,
    Indeterminate
  ];
