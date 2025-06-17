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
 * LAPACK based solver for all double complex eigenvalues/-vectors of a       *
 * hermitian matrix                                                           *
 *                                                                            *
 * -------------------------------------------------------------------------- */


#include "../inc.d/eigs.h"


// Eigenvalues and eigenvectors
void zheigsa(uint32_t n,
             const double complex *phi,
             bool evs,
             eigs_result *result) {

    // Copy matrix
    double complex *phi_cpy;
    phi_cpy = (double complex *)malloc(n*n*sizeof(double complex));
    for (uint32_t i=0; i<n*n; i++) phi_cpy[i] = phi[i];

    // Check if eigenvectors are desired
    char jobz;
    if (evs) {
        jobz = 'V';
    } else {
        jobz = 'N';
    }

    // Solve eigenproblem using LAPACK
    uint32_t i;
    double *eigvals = (double *)malloc(n*sizeof(double));
    lapack_int info = LAPACKE_zheev(LAPACK_ROW_MAJOR,
                                    jobz,
                                    'U',
                                    n,
                                    phi_cpy,
                                    n,
                                    eigvals);
    for (i=0; i<n; i++) result->eigvals[i] = CMPLX(eigvals[i], 0.);

    // Check result
    if (info) {
        printf("%s\n", "EIGS: LAPACKE_zheev did not converge"); exit(1);
    }

    // Check if eigenvectors are desired
    if (evs) {
        for (uint32_t i=0; i<n*n; i++) result->eigvecs[i] = phi_cpy[i];
    }

    // Clean up
    free(phi_cpy); free(eigvals);
}
