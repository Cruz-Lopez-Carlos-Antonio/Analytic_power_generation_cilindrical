(* ============================================================ *)
(* DIRECT AUGMENTED SOLVER                                      *)
(* Unknowns: {F(R), W(R), J(R), Omega}                          *)
(* ============================================================ *)

ClearAll[DirectSolver];

DirectSolver[nGrid_Integer?Positive] := Module[

    {
        PiD,
        K,
        h,
        rGrid,
        psiVals,
        sVals,
        mVals,

        idxF,
        idxW,
        idxJ,
        idxOmega,
        nUnknowns,

        entries,
        rhs,
        row,
        A,
        sol,

        Fvals,
        Wvals,
        Jvals,
        Pvals,
        OmegaDirect,

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
    (* Parameters                                               *)
    (* ======================================================== *)

    PiD = alpha*xi^2;
    K   = Lambda*delta^2;

    h = 1/nGrid;

    rGrid = N[Range[0, nGrid]*h];

    psiVals = PsiPB /@ rGrid;
    sVals   = Sinh /@ psiVals;
    mVals   = MR /@ rGrid;

    (* ======================================================== *)
    (* Unknown indexing                                         *)
    (* ======================================================== *)

    idxF[i_] := i + 1;
    idxW[i_] := nGrid + 2 + i;
    idxJ[i_] := 2*nGrid + 3 + i;

    idxOmega = 3*nGrid + 4;
    nUnknowns = 3*nGrid + 4;

    (* ======================================================== *)
    (* Matrix and RHS                                           *)
    (* ======================================================== *)

    entries = {};
    rhs = ConstantArray[0., nUnknowns];

    row = 0;

    (* ======================================================== *)
    (* F equations                                              *)
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
                 (rGrid[[i + 1]]/
                  mVals[[i + 1]])
        ];

        AppendTo[
            entries,
            {row, idxW[i + 1]} ->
                -(h/2)*
                 (rGrid[[i + 2]]/
                  mVals[[i + 2]])
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

    rhs[[row]] = 0;

    (* ======================================================== *)
    (* W equations                                              *)
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

    rhs[[row]] = PiD;

    (* ======================================================== *)
    (* J equations                                              *)
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

    rhs[[row]] = 0;

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

    rhs[[row]] = 0;

    (* ======================================================== *)
    (* Equation-count check                                     *)
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
    (* Solve global sparse system                               *)
    (* ======================================================== *)

    A = SparseArray[
        entries,
        {nUnknowns, nUnknowns}
    ];

    sol = LinearSolve[
        A,
        rhs
    ];

    (* ======================================================== *)
    (* Extract solution                                         *)
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

    (* Reconstruct P = R^2 W *)

    Pvals =
        rGrid^2*Wvals;

    (* ======================================================== *)
    (* Internal diagnostics                                     *)
    (* ======================================================== *)

    linearResidual =
        Max[
            Abs[
                A.sol - rhs
            ]
        ];

    (* Axis regularity *)

    axisResidual =
        Abs[
            2*Wvals[[1]]
            -
            PiD
            -
            OmegaDirect*K*sVals[[1]]
        ];

    (* F(1) = 0 *)

    FwallResidual =
        Abs[
            Fvals[[-1]]
        ];

    (* J(0) = 0 *)

    JaxisResidual =
        Abs[
            Jvals[[1]]
        ];

    (* J(1) = Omega *)

    selfConsistencyResidual =
        Abs[
            Jvals[[-1]]
            -
            OmegaDirect
        ];

    (* ======================================================== *)
    (* Discrete differential/integral residuals                 *)
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
        Max[Abs[residualF]];

    maxResidualW =
        Max[Abs[residualW]];

    maxResidualJ =
        Max[Abs[residualJ]];

    (* ======================================================== *)
    (* Complete radial table                                    *)
    (* ======================================================== *)

    solutionTable =
        Transpose[
            {
                rGrid,
                Fvals,
                Wvals,
                Pvals,
                Jvals
            }
        ];

    (* ======================================================== *)
    (* Return all information                                   *)
    (* ======================================================== *)

    <|

        "nGrid" -> nGrid,

        "PiD" -> PiD,

        "K" -> K,

        "Omega" -> OmegaDirect,

        "R" -> rGrid,

        "F" -> Fvals,

        "W" -> Wvals,

        "P" -> Pvals,

        "J" -> Jvals,

        "SolutionTable" -> solutionTable,

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