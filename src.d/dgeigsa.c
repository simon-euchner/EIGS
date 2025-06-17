/* -------------------------------------------------------------------------- *
 *                                                                            *
 * This file is part of the EIGS C-library by Simon Euchner.                  *
 *                                                                            *
 * -------------------------------------------------------------------------- *
 *                                                                            *
 * LICENSE: GPL-3.0                                                           *
 *                                                                            *
 * IMPORTANT: THIS IS FREE SOFTWARE WITHOUT ANY WARRANTY. THE USER IS FREE TO *
 *            MODIFY AND REDISTRIBUTE THIS SOFTWARE UNDER THE TERMS OF THE    *
 *            LICENSE LISTED ABOVE PUBLISHED BY THE FREE SOFTWARE FOUNDATION. *
 *            THE PUBLISHER, SIMON EUCHNER, IS NOT RESPONSIBLE FOR ANY        *
 *            NEGATIVE EFFECTS THIS SOFTWARE MAY CAUSE.                       *
 *                                                                            *
 * -------------------------------------------------------------------------- *
 *                                                                            *
 * LAPACK based solver for all double eigenvalues/-vectors                    *
 *                                                                            *
 * -------------------------------------------------------------------------- */


#include "../inc.d/eigs.h"


static void extract(eigs_result *, double *, double *, double *, bool);


// Eigenvalues and eigenvectors
void dgeigsa(uint32_t n,
             const double *phi,
             bool evs,
             eigs_result *result) {

    // Copy matrix
    double *phi_cpy;
    phi_cpy = (double *)malloc(n*n*sizeof(double));
    for (uint32_t i=0; i<n*n; i++) phi_cpy[i] = phi[i];


    // Solve eigenproblem using LAPACK
    double *wr, *wi, *vs;
    wr = (double *)malloc(n*sizeof(double));
    wi = (double *)malloc(n*sizeof(double));
    vs = (double *)malloc(n*n*sizeof(double));
    lapack_int info = LAPACKE_dgeev(LAPACK_ROW_MAJOR,
                                    'N',
                                    'V',
                                    n,
                                    phi_cpy,
                                    n,
                                    wr,
                                    wi,
                                    NULL,
                                    1,
                                    vs,
                                    n);

    // Check result
    if (info) {
        printf("%s\n", "EIGS: LAPACKE_dgeev did not converge"); exit(1);
    }

    // Extract eigenvalues and (possibly) eigenvectors
    extract(result, wr, wi, vs, evs);

    // Clean up
    free(wr); free(wi); free(vs); free(phi_cpy);
}


// Extract eigenvalues and (possiby) eigenvectors
static void extract(eigs_result *result,
                    double *wr,
                    double *wi,
                    double *vs,
                    bool evs) {

    uint32_t i, j, l, n = result->n;

    // Eigenvalues
    for (i=0; i<n; i++) result->eigvals[i] = CMPLX(wr[i], wi[i]);

    // Eigenvectors
    if (evs) {
        for (i=0; i<n; i++) {
            if (wi[i] == (double)0.) {
                for (j=0; j<n; j++)
                    result->eigvecs[j*n+i] = CMPLX(vs[j*n+i], 0.);
            } else {
                l = i++;
                for (j=0; j<n; j++) {
                    result->eigvecs[j*n+l] = CMPLX(vs[j*n+l],  vs[j*n+i]);
                    result->eigvecs[j*n+i] = CMPLX(vs[j*n+l], -vs[j*n+i]);
                }
            }
        }
    }
}
