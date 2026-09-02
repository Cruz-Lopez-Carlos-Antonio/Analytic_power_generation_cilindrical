(* ============================================================ *)
(* Code developed by Sánchez Lozano, G., Cruz-López C.-A., and  *)
(* F. Méndez, for the research:                                 *)
(* "Power generation in cylindrical microchannels with high     *)
(* surface zeta potential: a survey on viscoelectric effects    *)
(* over streaming potentials and efficiency"                    *)
(* Year: 2026                                                   *)
(* ------------------------------------------------------------ *)
(* DIRECT GENERAL MULTIPARAMETRIC SOLVER                        *)
(* Unknowns: {F(R), W(R), J(R), Omega}                          *)
(* Configurable inputs: PsiS, omega, Lambda, PiD                *)
(* ============================================================ *)

directSolverDirectory = DirectoryName[$InputFileName];

Scan[
    Get[
        FileNameJoin[
            {directSolverDirectory, #}
        ]
    ] &,
    {
        "Parameters.wl",
        "PoissonBoltzmann_parametric.wl",
        "F_cc_parametric.wl",
        "Lambda_parameter_parametric.wl",
        "MR_parametric.wl"
    }
];


ClearAll[DirectParametricSolver];


DirectParametricSolver[
    deltaInput_,
    config_ : <||>,
    nGrid_ : 1000
] := Module[

    {
        (* Parametric upstream quantities *)
        pb,
        cfg,
        psiSSetting,
        omegaSetting,
        lambdaSetting,
        piDSetting,
        deltaLoc,
        psiSLoc,
        omegaLoc,
        psiFun,
        centerPotential,
        FccValue,
        LambdaValue,

        (* Direct-solver coefficients *)
        PiD,
        K,
        h,
        rGrid,
        psiVals,
        sVals,
        mVals,

        (* Unknown indexing *)
        idxF,
        idxW,
        idxJ,
        idxOmega,
        nUnknowns,

        (* Matrix system *)
        entries,
        rhs,
        row,
        A,
        sol,

        (* Solution *)
        Fvals,
        Wvals,
        Jvals,
        Pvals,
        OmegaDirect,

        (* Phi_av quantities *)
        integralRF,
        GValue,
        dPhiDZValue,
        PhiAvValue,

        (* Diagnostics *)
        linearResidual,
        axisResidual,
        FwallResidual,
        JaxisResidual,
        selfConsistencyResidual,

        residualF,
        residualW,
        residualJ,

        maxResidualF,
        maxResidualW,
        maxResidualJ,

        solutionTable
    },


    (* ======================================================== *)
    (* 1. RESOLVE GENERAL CONFIGURATION                         *)
    (* ======================================================== *)

    cfg =
        Join[
            <|
                "PsiS"   -> Automatic,
                "omega"  -> Automatic,
                "Lambda" -> Automatic,
                "PiD"    -> Automatic
            |>,
            config
        ];

    psiSSetting   = cfg["PsiS"];
    omegaSetting  = cfg["omega"];
    lambdaSetting = cfg["Lambda"];
    piDSetting    = cfg["PiD"];


    psiSLoc =
        If[
            psiSSetting === Automatic,
            If[
                ValueQ[PsiS],
                PsiS,
                Print[
                    "ERROR: PsiS -> Automatic was requested, ",
                    "but the global symbol PsiS has no value."
                ];
                Return[$Failed]
            ],
            psiSSetting
        ];


    omegaLoc =
        If[
            omegaSetting === Automatic,
            If[
                ValueQ[omega],
                omega,
                Print[
                    "ERROR: omega -> Automatic was requested, ",
                    "but the global symbol omega has no value."
                ];
                Return[$Failed]
            ],
            omegaSetting
        ];


    If[
        !NumericQ[psiSLoc] || !NumericQ[omegaLoc],
        Print["ERROR: PsiS and omega must resolve to numeric values."];
        Return[$Failed];
    ];


    (* ======================================================== *)
    (* 2. PARAMETRIC UPSTREAM CHAIN                             *)
    (* ======================================================== *)

    pb =
        SolvePB[
            deltaInput,
            psiSLoc
        ];

    If[pb === $Failed,
        Print["ERROR: SolvePB failed."];
        Return[$Failed];
    ];


    deltaLoc = pb["Delta"];
    psiSLoc  = pb["PsiS"];

    psiFun =
        pb["Psi"];

    centerPotential =
        pb["CenterPotential"];


    FccValue =
        ComputeFcc[pb];

    If[!NumericQ[FccValue],
        Print["ERROR: ComputeFcc did not return a numeric value."];
        Return[$Failed];
    ];


    LambdaValue =
        If[
            lambdaSetting === Automatic,
            ComputeLambda[FccValue],
            lambdaSetting
        ];

    If[!NumericQ[LambdaValue],
        Print[
            "ERROR: Lambda did not resolve to a numeric value."
        ];
        Return[$Failed];
    ];


    (* ======================================================== *)
    (* 3. PARAMETERS OF THE DIRECT PROBLEM                      *)
    (* ======================================================== *)

    PiD =
        If[
            piDSetting === Automatic,
            alpha*xi^2,
            piDSetting
        ];

    If[!NumericQ[PiD],
        Print[
            "ERROR: PiD did not resolve to a numeric value."
        ];
        Return[$Failed];
    ];

    K =
        LambdaValue*deltaLoc^2;

    h =
        1/nGrid;

    rGrid =
        N[
            Range[0, nGrid]*h
        ];


    (* ======================================================== *)
    (* 4. NODAL COEFFICIENTS FOR CURRENT CONFIGURATION          *)
    (* ======================================================== *)

    psiVals =
        psiFun /@ rGrid;

    sVals =
        Sinh /@ psiVals;

    mVals =
        (
            ComputeM[
                pb,
                omegaLoc,
                #
            ]
        ) & /@ rGrid;


    If[
        !And @@ (NumericQ /@ mVals),

        Print[
            "ERROR: ComputeM did not return numeric values ",
            "on the complete radial grid."
        ];

        Return[$Failed];
    ];


    (* ======================================================== *)
    (* 4. UNKNOWN INDEXING                                      *)
    (* ======================================================== *)

    idxF[i_] :=
        i + 1;

    idxW[i_] :=
        nGrid + 2 + i;

    idxJ[i_] :=
        2*nGrid + 3 + i;

    idxOmega =
        3*nGrid + 4;

    nUnknowns =
        3*nGrid + 4;


    (* ======================================================== *)
    (* 5. MATRIX AND RHS                                        *)
    (* ======================================================== *)

    entries =
        {};

    rhs =
        ConstantArray[
            0.,
            nUnknowns
        ];

    row =
        0;


    (* ======================================================== *)
    (* 6. F EQUATIONS                                           *)
    (*                                                          *)
    (* F_{i+1} - F_i                                            *)
    (*   - h/2 [(R_i/M_i) W_i + (R_{i+1}/M_{i+1}) W_{i+1}] = 0 *)
    (* ======================================================== *)

    Do[

        row++;

        AppendTo[
            entries,
            {row, idxF[i]} -> -1
        ];

        AppendTo[
            entries,
            {row, idxF[i + 1]} -> 1
        ];

        AppendTo[
            entries,
            {row, idxW[i]} ->
                -(h/2)*
                 (
                    rGrid[[i + 1]]/
                    mVals[[i + 1]]
                 )
        ];

        AppendTo[
            entries,
            {row, idxW[i + 1]} ->
                -(h/2)*
                 (
                    rGrid[[i + 2]]/
                    mVals[[i + 2]]
                 )
        ];

        ,
        {i, 0, nGrid - 1}
    ];


    (* F(1) = 0 *)

    row++;

    AppendTo[
        entries,
        {row, idxF[nGrid]} -> 1
    ];

    rhs[[row]] =
        0;


    (* ======================================================== *)
    (* 7. W EQUATIONS                                           *)
    (*                                                          *)
    (* R_{i+1}^2 W_{i+1} - R_i^2 W_i                           *)
    (*  - h K/2 [R_i sinh(Psi_i)+R_{i+1}sinh(Psi_{i+1})]Omega  *)
    (*  = h PiD/2 (R_i + R_{i+1})                              *)
    (* ======================================================== *)

    Do[

        row++;

        AppendTo[
            entries,
            {row, idxW[i]} ->
                -rGrid[[i + 1]]^2
        ];

        AppendTo[
            entries,
            {row, idxW[i + 1]} ->
                rGrid[[i + 2]]^2
        ];

        AppendTo[
            entries,
            {row, idxOmega} ->
                -(h*K/2)*
                 (
                    rGrid[[i + 1]]*
                    sVals[[i + 1]]
                    +
                    rGrid[[i + 2]]*
                    sVals[[i + 2]]
                 )
        ];

        rhs[[row]] =
            (h*PiD/2)*
            (
                rGrid[[i + 1]]
                +
                rGrid[[i + 2]]
            );

        ,
        {i, 0, nGrid - 1}
    ];


    (* -------------------------------------------------------- *)
    (* Regular axis condition                                   *)
    (*                                                          *)
    (* 2 W(0) - Omega K sinh(Psi(0)) = PiD                     *)
    (* -------------------------------------------------------- *)

    row++;

    AppendTo[
        entries,
        {row, idxW[0]} -> 2
    ];

    AppendTo[
        entries,
        {row, idxOmega} ->
            -K*sVals[[1]]
    ];

    rhs[[row]] =
        PiD;


    (* ======================================================== *)
    (* 8. J EQUATIONS                                           *)
    (* ======================================================== *)

    Do[

        row++;

        AppendTo[
            entries,
            {row, idxJ[i]} -> -1
        ];

        AppendTo[
            entries,
            {row, idxJ[i + 1]} -> 1
        ];

        AppendTo[
            entries,
            {row, idxF[i]} ->
                -(h/2)*
                 rGrid[[i + 1]]*
                 sVals[[i + 1]]
        ];

        AppendTo[
            entries,
            {row, idxF[i + 1]} ->
                -(h/2)*
                 rGrid[[i + 2]]*
                 sVals[[i + 2]]
        ];

        ,
        {i, 0, nGrid - 1}
    ];


    (* J(0) = 0 *)

    row++;

    AppendTo[
        entries,
        {row, idxJ[0]} -> 1
    ];

    rhs[[row]] =
        0;


    (* J(1) - Omega = 0 *)

    row++;

    AppendTo[
        entries,
        {row, idxJ[nGrid]} -> 1
    ];

    AppendTo[
        entries,
        {row, idxOmega} -> -1
    ];

    rhs[[row]] =
        0;


    (* ======================================================== *)
    (* 9. EQUATION-COUNT CHECK                                  *)
    (* ======================================================== *)

    If[
        row =!= nUnknowns,

        Print[
            "ERROR: incorrect equation count. ",
            "Rows = ", row,
            ", unknowns = ", nUnknowns
        ];

        Return[$Failed];
    ];


    (* ======================================================== *)
    (* 10. SOLVE GLOBAL SPARSE SYSTEM                           *)
    (* ======================================================== *)

    A =
        SparseArray[
            entries,
            {nUnknowns, nUnknowns}
        ];

    sol =
        LinearSolve[
            A,
            rhs
        ];


    If[sol === $Failed,
        Print["ERROR: LinearSolve failed."];
        Return[$Failed];
    ];


    (* ======================================================== *)
    (* 11. EXTRACT SOLUTION                                     *)
    (* ======================================================== *)

    Fvals =
        sol[[
            idxF[0] ;; idxF[nGrid]
        ]];

    Wvals =
        sol[[
            idxW[0] ;; idxW[nGrid]
        ]];

    Jvals =
        sol[[
            idxJ[0] ;; idxJ[nGrid]
        ]];

    OmegaDirect =
        sol[[idxOmega]];


    Pvals =
        rGrid^2*Wvals;


    (* ======================================================== *)
    (* 12. COMPUTE G, dPhi/dZ AND Phi_av                        *)
    (* ======================================================== *)

    integralRF =
        Sum[
            (h/2)*
            (
                rGrid[[i + 1]]*
                Fvals[[i + 1]]
                +
                rGrid[[i + 2]]*
                Fvals[[i + 2]]
            ),
            {i, 0, nGrid - 1}
        ];


    If[
        !NumericQ[integralRF] || PossibleZeroQ[integralRF],

        Print[
            "ERROR: the numerical integral Integral_0^1 R F(R)dR ",
            "is zero or nonnumeric."
        ];

        Return[$Failed];
    ];


    GValue =
        1/(2*integralRF);


    dPhiDZValue =
        -LambdaValue*
        deltaLoc^2*
        GValue*
        OmegaDirect;


    PhiAvValue =
        (1/2)*
        dPhiDZValue;


    (* ======================================================== *)
    (* 13. INTERNAL DIAGNOSTICS                                 *)
    (* ======================================================== *)

    linearResidual =
        Max[
            Abs[
                A.sol - rhs
            ]
        ];


    axisResidual =
        Abs[
            2*Wvals[[1]]
            -
            PiD
            -
            OmegaDirect*K*sVals[[1]]
        ];


    FwallResidual =
        Abs[
            Fvals[[-1]]
        ];


    JaxisResidual =
        Abs[
            Jvals[[1]]
        ];


    selfConsistencyResidual =
        Abs[
            Jvals[[-1]]
            -
            OmegaDirect
        ];


    (* ======================================================== *)
    (* 14. DISCRETE EQUATION RESIDUALS                          *)
    (* ======================================================== *)

    residualF =
        Table[
            Fvals[[i + 2]]
            -
            Fvals[[i + 1]]
            -
            (h/2)*
            (
                rGrid[[i + 1]]*
                Wvals[[i + 1]]/
                mVals[[i + 1]]
                +
                rGrid[[i + 2]]*
                Wvals[[i + 2]]/
                mVals[[i + 2]]
            ),
            {i, 0, nGrid - 1}
        ];


    residualW =
        Table[
            rGrid[[i + 2]]^2*
            Wvals[[i + 2]]
            -
            rGrid[[i + 1]]^2*
            Wvals[[i + 1]]
            -
            (h*K/2)*
            (
                rGrid[[i + 1]]*
                sVals[[i + 1]]
                +
                rGrid[[i + 2]]*
                sVals[[i + 2]]
            )*
            OmegaDirect
            -
            (h*PiD/2)*
            (
                rGrid[[i + 1]]
                +
                rGrid[[i + 2]]
            ),
            {i, 0, nGrid - 1}
        ];


    residualJ =
        Table[
            Jvals[[i + 2]]
            -
            Jvals[[i + 1]]
            -
            (h/2)*
            (
                rGrid[[i + 1]]*
                sVals[[i + 1]]*
                Fvals[[i + 1]]
                +
                rGrid[[i + 2]]*
                sVals[[i + 2]]*
                Fvals[[i + 2]]
            ),
            {i, 0, nGrid - 1}
        ];


    maxResidualF =
        Max[
            Abs[residualF]
        ];

    maxResidualW =
        Max[
            Abs[residualW]
        ];

    maxResidualJ =
        Max[
            Abs[residualJ]
        ];


    (* ======================================================== *)
    (* 15. COMPLETE RADIAL TABLE                                *)
    (* ======================================================== *)

    solutionTable =
        Transpose[
            {
                rGrid,
                psiVals,
                mVals,
                Fvals,
                Wvals,
                Pvals,
                Jvals
            }
        ];


    (* ======================================================== *)
    (* 16. RETURN COMPLETE PARAMETRIC SOLUTION                  *)
    (* ======================================================== *)

    <|

        "Delta" ->
            deltaLoc,

        "Config" ->
            cfg,

        "PsiS" ->
            psiSLoc,

        "omega" ->
            omegaLoc,

        "CenterPotential" ->
            centerPotential,

        "nGrid" ->
            nGrid,

        "Fcc" ->
            FccValue,

        "Lambda" ->
            LambdaValue,

        "LambdaMode" ->
            If[lambdaSetting === Automatic, "Automatic", "Imposed"],

        "PiD" ->
            PiD,

        "PiDMode" ->
            If[piDSetting === Automatic, "Automatic", "Imposed"],

        "K" ->
            K,

        "Omega" ->
            OmegaDirect,

        "IntegralRF" ->
            integralRF,

        "G" ->
            GValue,

        "dPhiDZ" ->
            dPhiDZValue,

        "PhiAv" ->
            PhiAvValue,

        "R" ->
            rGrid,

        "Psi" ->
            psiVals,

        "M" ->
            mVals,

        "F" ->
            Fvals,

        "W" ->
            Wvals,

        "P" ->
            Pvals,

        "J" ->
            Jvals,

        "SolutionTable" ->
            solutionTable,

        "LinearResidual" ->
            linearResidual,

        "AxisResidual" ->
            axisResidual,

        "FwallResidual" ->
            FwallResidual,

        "JaxisResidual" ->
            JaxisResidual,

        "SelfConsistencyResidual" ->
            selfConsistencyResidual,

        "MaxResidualF" ->
            maxResidualF,

        "MaxResidualW" ->
            maxResidualW,

        "MaxResidualJ" ->
            maxResidualJ

    |>

];


(* ============================================================ *)