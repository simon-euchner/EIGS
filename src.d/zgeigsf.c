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
 * ARPACK based solver for a few eigenvalues/-vectors of a general double     *
 * complex endomorphism                                                       *
 * -------------------------------------------------------------------------- */


#include "../inc.d/eigs.h"


// Data for internal usage
typedef struct _ZgeigsfData {

    // User set
    a_int n;
    zeigs_phi *phi;
    void *phi_data;
    a_int nev;
    const char *which;
    bool evs;
    double tol;
    a_int ncv;
    a_int mxiter;

    // Internal
    a_int ido;
    const char *bmat;
    a_dcomplex *resid;
    a_dcomplex *v;
    a_int ldv;
    a_int *iparam;
    a_int *ipntr;
    a_dcomplex *workd;
    a_int lworkl;
    a_dcomplex *workl;
    double *rwork;
    a_int info;
    a_int ldz;
    a_dcomplex *workev;

    // Results
    a_dcomplex *d;
    a_dcomplex *z;

} zgeigsf_data;


static zgeigsf_data *zgeigsf_init(a_int,
                                  zeigs_phi *,
                                  void *,
                                  a_int,
                                  const char *,
                                  bool,
                                  double,
                                  a_int);
static void zgeigsf_data_destroy(zgeigsf_data *);
static void arnoldi_iterations(zgeigsf_data *);
static void iterate(zgeigsf_data *);
static void extract(zgeigsf_data *);
static eigs_result *prepare_result(zgeigsf_data *, eigs_result *);


// Eigenvalues and eigenvectors
void zgeigsf(a_int n,
             zeigs_phi *phi,
             void *phi_data,
             bool evs,
             const char *which,
             a_int k,
             double tol,
             a_int maxiter,
             eigs_result *result) {

    // Initialize data
    zgeigsf_data *data = zgeigsf_init(n,
                                      phi,
                                      phi_data,
                                      k,
                                      which,
                                      evs,
                                      tol,
                                      maxiter);

    // Arnoldi iterations
    arnoldi_iterations(data);

    // Extract eigenvalues and (possibly) eigenvectors
    extract(data);

    // Prepare result
    prepare_result(data, result);

    // Clean up
    zgeigsf_data_destroy(data);
}

// Initialize eigenproblem
static zgeigsf_data *zgeigsf_init(a_int n,
                                  zeigs_phi *phi,
                                  void *phi_data,
                                  a_int k,
                                  const char *which,
                                  bool evs,
                                  double tol,
                                  a_int maxiter) {

    // Allocate memory for data
    zgeigsf_data *data = (zgeigsf_data *)malloc(sizeof(zgeigsf_data));

    // User set
    data->n = n;
    data->phi = phi;
    data->nev = k;
    data->which = which;
    data->evs = evs;
    data->tol = tol; // Default 0. (machine precision)
    data->mxiter = maxiter; // Default 10*n
    data->phi_data = phi_data; // Default NULL

    // Internal
    data->ido = 0;
    data->bmat = "I";
    data->resid = (a_dcomplex *)malloc(n*sizeof(a_dcomplex));
    if ((data->ncv = 2*k+1) < 20) data->ncv = 20;
    if (data->ncv > n) data->ncv = n;
    data->v = (a_dcomplex *)calloc(n*data->ncv, sizeof(a_dcomplex));
    data->ldv = n;
    data->iparam = (a_int *)calloc(11, sizeof(a_int));
    data->iparam[0] = 1;
    data->iparam[2] = maxiter;
    data->iparam[3] = 1;
    data->iparam[6] = 1;
    data->ipntr = (a_int *)calloc(14, sizeof(a_int));
    data->workd = (a_dcomplex *)calloc(3*n, sizeof(a_dcomplex));
    data->lworkl = 3*data->ncv*(data->ncv+2);
    data->workl = (a_dcomplex *)calloc(data->lworkl, sizeof(a_dcomplex));
    data->rwork = (double *)calloc(data->ncv, sizeof(double));
    data->info = 0;
    data->ldz = n;
    data->workev = (a_dcomplex *)calloc(3*data->ncv, sizeof(a_dcomplex));

    // Results
    data->d = (a_dcomplex *)calloc(data->nev+1, sizeof(a_dcomplex));
    data->z = (a_dcomplex *)calloc(n*data->nev, sizeof(a_dcomplex));

    return data;
}

// Free for zeigsf_data type
static void zgeigsf_data_destroy(zgeigsf_data *data) {
    free(data->resid); data->resid = NULL;
    free(data->v); data->v = NULL;
    free(data->iparam); data->iparam = NULL;
    free(data->ipntr); data->ipntr = NULL;
    free(data->workd); data->workd = NULL;
    free(data->workl); data->workl = NULL;
    free(data->rwork); data->rwork = NULL;
    free(data->workev); data->workev = NULL;
    free(data->d); data->d = NULL;
    free(data->z); data->z = NULL;
    free(data);
}

// Do Arnoldi iterations
static void arnoldi_iterations(zgeigsf_data *data) {

    // Arnoldi iterations
    do {
        iterate(data);
    } while ((data->ido == 1) || (data->ido == -1));

    // Check for errors
    if (data->ido != 99) {
        printf("%s\n", "ZEIGSF: ARNOLDI PROCESS DID NOT CONVERGE");
        exit(1);
    }
}

// Do a single Arnoldi iteration
static void iterate(zgeigsf_data *data) {

    // Call ZNAUPD
    znaupd_c(&data->ido,
             data->bmat,
             data->n,
             data->which,
             data->nev,
             data->tol,
             data->resid,
             data->ncv,
             data->v,
             data->ldv,
             data->iparam,
             data->ipntr,
             data->workd,
             data->workl,
             data->lworkl,
             data->rwork,
             &data->info);

    // Check for errors
    int nerror = 0;
    if ((data->ido != 1) && (data->ido != -1) && (data->ido != 99)) {
        printf("ZEIGSF: ERROR DURING ITERATION: IDO = %d\n", data->ido);
        nerror++;
    }
    if ((data->info != 0) && (data->info != 1)) {
        printf("ZEIGSF: ERROR DURING ITERATION: INFO = %d\n", data->info);
        nerror++;
    }
    if (data->info == 1) {
        printf("%s\n", "ZEIGSF: MAXIMAL ALLOWED ITERATIONS REACHED");
        nerror++;
    }
    if (nerror) exit(1);

    // Compute action of phi
    a_int xpntr = data->ipntr[0]-1;
    a_int ypntr = data->ipntr[1]-1;
    data->phi(data->phi_data,
              data->n,
              &(data->workd[xpntr]),
              &(data->workd[ypntr]));
}

// Extract eigenvalues and (possiby) eigenvectors
static void extract(zgeigsf_data *data) {

    // For internal use
    const char *howmny = "A"; if (data->evs) howmny = "P";
    a_int *select = (a_int *)calloc(data->ncv, sizeof(a_int));
    a_dcomplex sigma = CMPLX(0., 0.); // Not referenced

    // Call ZNEUPD
    zneupd_c(data->evs,
             howmny,
             select,
             data->d,
             data->z,
             data->ldz,
             sigma,
             data->workev,
             data->bmat,
             data->n,
             data->which,
             data->nev,
             data->tol,
             data->resid,
             data->ncv,
             data->v,
             data->ldv,
             data->iparam,
             data->ipntr,
             data->workd,
             data->workl,
             data->lworkl,
             data->rwork,
             &data->info);

    // Clean up
    free(select);

    // Check for errors
    if (data->info) {
        printf("ZEIGSF: COULD NOT EXTRACT RESULTS: INFO = %d\n", data->info);
        exit(1);
    }
}

// Load data into result and reorder it to row major
static eigs_result *prepare_result(zgeigsf_data *data, eigs_result *result) {
    a_int n, k, i, j, count;
    n = data->n; k = data->nev; count = 0;
    result->n = n; result->k = k;
    for (j=0; j<k; j++) result->eigvals[j] = data->d[j];
    if (data->evs) {
        for (i=0; i<n; i++) {
            for (j=0; j<k; j++) {
                result->eigvecs[count++] = data->z[n*j+i];
            }
        }
    } else {
        result->eigvecs = NULL;
    }
    return result;
}
