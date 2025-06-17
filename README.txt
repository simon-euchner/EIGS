o ---------------------------------------------------------------------------- o
| Eigensolver for the C programming language by Simon Euchner                  |
o ---------------------------------------------------------------------------- o


LICENSE NOTICE.

    LICENSES: GPL-3.0

    IMPORTANT: THIS IS FREE SOFTWARE WITHOUT ANY WARRANTY. THE USER IS FREE TO
               MODIFY AND REDISTRIBUTE THIS SOFTWARE UNDER THE TERMS OF THE
               LICENSE LISTED ABOVE PUBLISHED BY THE FREE SOFTWARE FOUNDATION.
               THE PUBLISHER, SIMON EUCHNER, IS NOT RESPONSIBLE FOR ANY NEGATIVE
               EFFECTS THIS SOFTWARE MAY CAUSE.

               LICENCES FOR EXTERNAL SOFTWARE ARE PROVIDED IN THE CORRESPONDING
               DIRECTORY(IES).

               EXTERNAL LINKS.

               SCIPY     : [1]
               ARPACK-NG : [2]


Prerequisites.

    - FORTRAN and C compiler

    - BLAS

    - LAPACK

    - LAPACKE


Functionality.

    There is only a single function, named "eigs", which is suited for any type
    of standard eigenproblem. The function "eigs" has the following sinature:

    eigs_result *eigs( const char                *solver            ,
                       zeigs_phi                 *zphi              ,
                       deigs_phi                 *dphi              ,
                       const double complex      *zphi_matrix       ,
                       const double              *dphi_matrix       ,
                       void                      *phi_data          ,
                       int32_t                    n                 ,
                       int32_t                    k                 ,
                       const char                *which             ,
                       int32_t                    maxiter           ,
                       double                     tol               ,
                       bool                       evs                 );

    --- Arguments. ---

    "solver": Type of input matrix.
                  "zg": general double complex
                  "dg": general double
                  "zh": hermitian double complex
                  "ds": symmetric double

    "zphi": Linear map to be diagonalized (Use if only a few eigenvalues/
            -vectors are desired).
                In this case "solver" must be either "zg" or "zh".
                For details on the type "zeigs_phi" see "./inc.d/eigs.h".
                If not used, pass "NULL".

    "dphi": Linear map to be diagonalized (Use if only a few eigenvalues/
            -vectors are desired).
                In this case "solver" must be either "dg" or "ds".
                For details on the type "deigs_phi" see "./inc.d/eigs.h".
                If not used, pass "NULL".

    "zphi_matrix": Row-major ordered double complex array that represents the
                   matrix to be diagoanlized (Use if all eigenvalues/-vectors
                   are desired).
                       In this case "solver" must be either "zg" or "zh".
                       If not used, pass "NULL".

    "dphi_matrix": Row-maojr ordered double array that represents the matrix to
                   be diagoanlized (Use if all eigenvalues/-vectors are
                   desired).
                       In this case "solver" must be either "dg" or "ds".
                       If not used, pass "NULL".

    "phi_data": Additional data that is passed to "zphi"("dphi").
                For details on the type "zeigs_phi"("deigs_phi") see
                "./inc.d/eigs.h".
                If not used, pass "NULL".

    "n": Dimension of the vector space on which the eigenproblem is formulated.

    "k": Number of desired eigenvalues/-vectors.
             If one either "zphi" or "dphi" is not "NULL", "k" must be between
             (including boundaries) 1 and n-2. If one of "zphi_matrix" and
             "dphi_matrix" is not "NULL", "k" must be equal to "n".

    "which": Which eigenvalues/-vectors to compute.
                 Only applies if a few egenvalues are desired. May be set to
                 "NULL" otherwise. Options are dependent on the type of input
                 matrix. See source code in "./ARPACK/SRC/" [2]. The flags "SM"
                 (smallest magnitude) and "LM" (largest magniude) are always
                 available.

    "maxiter": Maximal number of allowed Arnoldi iterations.
                   Only applies of either "zphi" or "dphi" is not "NULL".
                   Use "-1" for defualt value (see also [1]).

    "tol": Relative tolerance of Ritz values.
               Only applies if either "zphi" or "dphi" is not "NULL".
               Use "-1." for default value "0.0".
               More information can be found in the source code in
               "./ARPACK/SRC/" [2].

    "evs": Decides if the eigenvectors are computed/returned or not.

    --- Return. ---

        The result is of type "eigs_result", a.k.a. "_EigsResult".
        This structure hosts four members:
            "n": Dimension of the vector space on which the eigenproblem is
                 formulated.
            "k": Number of desired eigenvalues/-vectors.
            "eigvals": Double complex array containing the eigenvalues.
            "eigvecs": Double complex array containing the eigenvectors, where
                       "(eigvecs[i*n+j], i=1,...,n)" is the j-th of the in total
                       "k" eigenvectors.

    --- Additional information. ---

        Only one of "zphi", "dphi", "zphi_matrix", and "dphi_matrix" can be not
        equal to "NULL".


General information.

    To keep things simple, I chose to always return the eigenvalues and
    eigenvectors in double precision complex numbers. If your problem is real,
    rest assured that, internally, "eigs" does not use the double complex
    routines. The type of solver is selected via the argument "solver" given to
    "eigs".

    THIS IS WORK IN PROGRESS, THE "REAL" SOLVER DO NOT WORK RIGHT YET !!!


Installation.

    Run "./Makefile" to generate the library "libeigs.so" at "./lib.d/". Place
    this library as well as the header "./inc.d/eigs.h" at a preferred location.
    Link with "-leigs -larpack -llapack -llapacke -lblas" and (possibly) add
    runtime dependencies using "-Wl,-rpath,<put-the-paths-here>". Also, do not
    forgest to tell the compiler where the libraries are located, using
    (possibly multiple) flags like "-L<a-path>".


External links.

    [1] https://www.gitub.com/scipy/scipy

    [2] https://www.gitub.com/opencollab/arpack-ng
