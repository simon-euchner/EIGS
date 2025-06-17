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
 * LAPACK based solver for all double complex eigenvalues/-vectors            *
 *                                                                            *
 * -------------------------------------------------------------------------- */


#include "../inc.d/eigs.h"


// Eigenvalues and eigenvectors
void zgeigsa(uint32_t n,
             const double complex *phi,
             bool evs,
             eigs_result *result) {

    // Copy matrix
    double complex *phi_cpy, *eigvecs;
    phi_cpy = (double complex *)malloc(n*n*sizeof(double complex));
    for (uint32_t i=0; i<n*n; i++) phi_cpy[i] = phi[i];

    // Check if eigenvectors are desired
    if (evs) {
        eigvecs = result->eigvecs;
    } else {
        eigvecs = (double complex *)malloc(n*n*sizeof(double complex));
    }


    // Solve eigenproblem using LAPACK
    lapack_int info = LAPACKE_zgeev(LAPACK_ROW_MAJOR,
                                    'N',
                                    'V',
                                    n,
                                    phi_cpy,
                                    n,
                                    result->eigvals,
                                    NULL,
                                    1,
                                    eigvecs,
                                    n);

    // Check result
    if (info) {
        printf("%s\n", "EIGS: LAPACKE_zgeev did not converge"); exit(1);
    }

    // Clean up
    free(phi_cpy); if (!evs) free(eigvecs);
}
