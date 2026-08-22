(* ============================================================ *)
(* G_parametric.wl                                              *)
(* Parametric implementation of G                               *)
(* ============================================================ *)
(*
   Computes

       IF = Integral_0^1 R F(R) dR

   and

       G = 1/(2 IF)

   from an already constructed parametric F object.

   The numerical strategy mirrors the validated G.wl module.

   IMPORTANT:
   No artificial SetPrecision is introduced.
*)
(* ============================================================ *)


ClearAll[ComputeGParametric];


ComputeGParametric[
   fObj_Association
   ] :=
 Module[
  {
   flowIntegrandParametric,
   IFParametric,
   GParametric
   },


  (* ========================================================= *)
  (* NUMERICAL SETTINGS                                       *)
  (* Same strategy as validated G.wl                           *)
  (* ========================================================= *)

  wpGParametric = 30;
  agGParametric = 12;
  pgGParametric = 12;
  maxRecursionGParametric = 30;


  (* ========================================================= *)
  (* FLOW INTEGRAND                                           *)
  (* ========================================================= *)

  ClearAll[flowIntegrandParametric];

  flowIntegrandParametric[R_?NumericQ] :=
   R*
    fObj["F"][R];


  (* ========================================================= *)
  (* FLOW INTEGRAL                                            *)
  (* ========================================================= *)

  IFParametric =
   NIntegrate[
    flowIntegrandParametric[R],

    {R, 0, 1},

    WorkingPrecision -> wpGParametric,
    AccuracyGoal -> agGParametric,
    PrecisionGoal -> pgGParametric,
    MaxRecursion -> maxRecursionGParametric,

    Method -> {
      "GlobalAdaptive",
      "SymbolicProcessing" -> 0
      }
    ];


  (* ========================================================= *)
  (* NON-DEGENERACY CHECK                                     *)
  (* ========================================================= *)

  If[
   Abs[IFParametric] < 10^-12,

   Print[
    "WARNING in G_parametric.wl: ",
    "Integral R F(R) dR is numerically small: ",
    N[IFParametric, 18]
    ]
   ];


  (* ========================================================= *)
  (* PRESSURE GRADIENT                                        *)
  (* ========================================================= *)

  GParametric =
   1/
    (
     2 IFParametric
     );


  (* ========================================================= *)
  (* RETURN OBJECT                                            *)
  (* ========================================================= *)

  <|
   "IF" -> IFParametric,
   "G" -> GParametric
   |>
  ];