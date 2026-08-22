(* ========================================================= *)
(* POISSON-BOLTZMANN PARAMETRIC SOLVER                       *)
(* ========================================================= *)
(*
   Solves the cylindrical nonlinear Poisson-Boltzmann
   equation

       Psi''(R) + (1/R) Psi'(R)
           = delta^2 Sinh[Psi(R)]

   for externally supplied values of

       deltaInput
       psiSInput

   with boundary conditions

       Psi'(0) = 0
       Psi(1)  = psiSInput.

   The solver is executed ONCE for each pair

       (deltaInput, psiSInput),

   and returns an Association containing the resulting
   functions Psi(R), Psi'(R), and p(R)=R Psi'(R).

   This avoids solving the shooting problem repeatedly
   every time Psi is evaluated downstream.
*)
(* ========================================================= *)


(* ========================================================= *)
(* 1. NUMERICAL SETTINGS                                     *)
(* ========================================================= *)

wpPBParam = 30;
wpExtPBParam = 40;

epsPBParam = SetPrecision[10^-8, wpExtPBParam];


(* ========================================================= *)
(* 2. PARAMETRIC SOLVER                                      *)
(* ========================================================= *)

ClearAll[SolvePB];

SolvePB[
   deltaInput_?NumericQ,
   psiSInput_?NumericQ
   ] :=
 Module[
  {
   deltaLoc,
   psiSLoc,
   psi,
   p,
   r,
   pbForCenter,
   wallResidual,
   psi0Guess,
   psi0,
   centerPotential,
   pbSolution,
   psiFunction,
   pFunction,
   psiPB,
   pPB,
   psiPrimePB
   },


  (* ------------------------------------------------------- *)
  (* PARAMETER PRECISION                                     *)
  (* ------------------------------------------------------- *)

  deltaLoc =
   SetPrecision[deltaInput, wpExtPBParam];

  psiSLoc =
   SetPrecision[psiSInput, wpExtPBParam];


  (* ------------------------------------------------------- *)
  (* INTEGRATION FOR A GIVEN CENTER POTENTIAL                *)
  (* ------------------------------------------------------- *)

  pbForCenter[psi0_?NumericQ] :=
   Module[
    {psi0Loc},

    psi0Loc =
     SetPrecision[psi0, wpExtPBParam];

    NDSolveValue[
     {
      psi'[r] == p[r]/r,

      p'[r] ==
       deltaLoc^2*r*Sinh[psi[r]],

      psi[epsPBParam] ==
       psi0Loc
        +
        (deltaLoc^2*Sinh[psi0Loc]/4)*
         epsPBParam^2,

      p[epsPBParam] ==
       (deltaLoc^2*Sinh[psi0Loc]/2)*
        epsPBParam^2
      },

     {psi, p},

     {r, epsPBParam, 1},

     WorkingPrecision -> wpPBParam,
     AccuracyGoal -> 18,
     PrecisionGoal -> 18,
     MaxSteps -> Infinity
     ]
    ];


  (* ------------------------------------------------------- *)
  (* SHOOTING RESIDUAL                                       *)
  (* ------------------------------------------------------- *)

  wallResidual[psi0_?NumericQ] :=
   Module[
    {sol},

    sol = pbForCenter[psi0];

    sol[[1]][1] - psiSLoc
    ];


  (* ------------------------------------------------------- *)
  (* DEBYE-HUCKEL INITIAL ESTIMATE                           *)
  (* ------------------------------------------------------- *)

  psi0Guess =
   psiSLoc/BesselI[0, deltaLoc];


  (* ------------------------------------------------------- *)
  (* CENTER POTENTIAL                                        *)
  (* ------------------------------------------------------- *)

  centerPotential =
   psi0 /. FindRoot[
      wallResidual[psi0] == 0,

      {
       psi0,
       (9/10)*psi0Guess,
       (11/10)*psi0Guess
       },

      WorkingPrecision -> wpPBParam,
      AccuracyGoal -> 20,
      PrecisionGoal -> 20,
      MaxIterations -> 100
      ];


  (* ------------------------------------------------------- *)
  (* FINAL NUMERICAL SOLUTION                                *)
  (* ------------------------------------------------------- *)

  pbSolution =
   pbForCenter[centerPotential];

  psiFunction =
   pbSolution[[1]];

  pFunction =
   pbSolution[[2]];


  (* ------------------------------------------------------- *)
  (* FUNCTIONS INCLUDING R = 0                               *)
  (* ------------------------------------------------------- *)

  psiPB[rIn_?NumericQ] :=
   Which[

    rIn == 0,
    centerPotential,

    0 < rIn < epsPBParam,
    centerPotential
     +
     (deltaLoc^2*Sinh[centerPotential]/4)*
      rIn^2,

    epsPBParam <= rIn <= 1,
    psiFunction[rIn],

    True,
    Indeterminate
    ];


  pPB[rIn_?NumericQ] :=
   Which[

    rIn == 0,
    0,

    0 < rIn < epsPBParam,
    (deltaLoc^2*Sinh[centerPotential]/2)*
     rIn^2,

    epsPBParam <= rIn <= 1,
    pFunction[rIn],

    True,
    Indeterminate
    ];


  psiPrimePB[rIn_?NumericQ] :=
   Which[

    rIn == 0,
    0,

    0 < rIn < epsPBParam,
    (deltaLoc^2*Sinh[centerPotential]/2)*
     rIn,

    epsPBParam <= rIn <= 1,
    pFunction[rIn]/rIn,

    True,
    Indeterminate
    ];


  (* ------------------------------------------------------- *)
  (* RETURN COMPLETE SOLUTION OBJECT                         *)
  (* ------------------------------------------------------- *)

  <|
   "Delta" -> deltaLoc,
   "PsiS" -> psiSLoc,
   "CenterPotential" -> centerPotential,
   "Psi" -> psiPB,
   "PsiPrime" -> psiPrimePB,
   "p" -> pPB
   |>
  ];