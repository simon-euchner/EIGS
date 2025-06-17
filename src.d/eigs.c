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
 * Main program: Implementation of the eigensolver *eigs*.                    *
 *                                                                            *
 * -------------------------------------------------------------------------- */


#include "../inc.d/eigs.h"


static eigs_result *eigs_result_alloc(int32_t, int32_t, bool);


// Eigensolver
eigs_result *eigs(const char *solver,
                  zeigs_phi *zphi,
                  deigs_phi *dphi,
                  const double complex *zphi_matrix,
                  const double *dphi_matrix,
                  void *phi_data,
                  int32_t n,
                  int32_t k,
                  const char *which,
                  int32_t maxiter,
                  double tol,
                  bool evs) {

    // Allocate memory for result
    eigs_result *result = eigs_result_alloc(n, k, evs);

    // Apply solver to problem
    if (!strcmp(solver, "zg")) { /* --- DOUBLE COMPLEX GENERAL --- */

        // Apply defaults if nessesary
        if (tol < 0.) tol = 0.;
        if (maxiter <= 0) maxiter = 10*n;

        // Either solve for all or a few eigenvalues/-vectors
        if (k == n) {
            // LAPACK(E)'s (LAPACK_)ZGEEV
            (void)zphi; (void)dphi; (void)dphi_matrix; (void)which;
            (void)maxiter; (void)tol; (void)evs;
            zgeigsa(n, zphi_matrix, evs, result);
        } else {
            // ARPACK's ZNAUPD and ZNEUPD (Carefull, make sure k < n-1!)
            (void)dphi; (void)zphi_matrix; (void)dphi_matrix;
            zgeigsf(n, zphi, phi_data, evs, which, k, tol, maxiter, result);
        }

    } else
    if (!strcmp(solver, "dg")) { /* --- DOUBLE GENERAL --- */

        // Apply defaults if nessesary
        if (tol < 0.) tol = 0.;
        if (maxiter <= 0) maxiter = 10*n;

        // Either solve for all or a few eigenvalues/-vectors
        if (k == n) {
            // LAPACK(E)'s (LAPACK_)DGEEV
            (void)zphi; (void)dphi; (void)zphi_matrix; (void)which;
            (void)maxiter; (void)tol; (void)evs;
            dgeigsa(n, dphi_matrix, evs, result);
        } else {
            // ARPACK's DNAUPD and DNEUPD (Carefull, make sure k < n-1!)
            (void)dphi; (void)zphi_matrix; (void)dphi_matrix;
            dgeigsf(n, dphi, phi_data, evs, which, k, tol, maxiter, result);
        }

    } else
    if (!strcmp(solver, "zh")) { /* --- DOUBLE COMPLEX HERMITIAN --- */

        // Apply defaults if nessesary
        if (tol < 0.) tol = 0.;
        if (maxiter <= 0) maxiter = 10*n;

        // Either solve for all or a few eigenvalues/-vectors
        if (k == n) {
            // LAPACK(E)'s (LAPACK_)DGEEV
            (void)zphi; (void)dphi; (void)dphi_matrix; (void)which;
            (void)maxiter; (void)tol; (void)evs;
            zheigsa(n, zphi_matrix, evs, result);
        } else {
            // ARPACK's ZNAUPD and ZNEUPD (Carefull, make sure k < n-1!)
            (void)zphi_matrix; (void)dphi; (void)dphi_matrix;
            zgeigsf(n, zphi, phi_data, evs, which, k, tol, maxiter, result);
        }

    } else
    if (!strcmp(solver, "ds")) { /* --- DOUBLE SYMMETRIC --- */

        // Apply defaults if nessesary
        if (tol < 0.) tol = 0.;
        if (maxiter <= 0) maxiter = 10*n;

        // Either solve for all or a few eigenvalues/-vectors
        if (k == n) {
            // LAPACK(E)'s (LAPACK_)DGEEV
            (void)zphi; (void)dphi; (void)zphi_matrix; (void)which;
            (void)maxiter; (void)tol; (void)evs;
            dseigsa(n, dphi_matrix, evs, result);
        } else {
            // ARPACK's ?NAUPD and ?NEUPD (Carefull, make sure k < n-1!)
            printf("%s\n", "Not implemented yet");
        }

    } else {

        printf("EIGS: Solver *%s* not implemented\n", solver);
        exit(1);

    }

    return result;
}

// Allocater for result type
static eigs_result *eigs_result_alloc(int32_t n, int32_t k, bool evs) {
    eigs_result *result = (eigs_result *)malloc(sizeof(eigs_result));
    result->n = n; result->k = k;
    result->eigvals = (double complex *)malloc(k*sizeof(double complex));
    if (evs)
        result->eigvecs = (double complex *)malloc(n*k*sizeof(double complex));
    else
        result->eigvecs = NULL;
    return result;
}

// Free memory allocated by result
void eigs_result_free(eigs_result *result) {
    free(result->eigvals);
    if (result->eigvecs) free(result->eigvecs);
    free(result);
}
