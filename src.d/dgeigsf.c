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
 * endomorphism                                                               *
 * -------------------------------------------------------------------------- */


#include "../inc.d/eigs.h"


// Data for internal usage
typedef struct _DgeigsfData {

    // User set
    a_int n;
    deigs_phi *phi;
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
    double *resid;
    double *v;
    a_int ldv;
    a_int *iparam;
    a_int *ipntr;
    double *workd;
    a_int lworkl;
    double *workl;
    a_int info;
    a_int ldz;
    double *workev;

    // Results
    double *dr;
    double *di;
    double *z;

} dgeigsf_data;


static dgeigsf_data *dgeigsf_init(a_int,
                                  deigs_phi *,
                                  void *,
                                  a_int,
                                  const char *,
                                  bool,
                                  double,
                                  a_int);
static void dgeigsf_data_destroy(dgeigsf_data *);
static void arnoldi_iterations(dgeigsf_data *);
static void iterate(dgeigsf_data *);
static void extract(dgeigsf_data *);
static eigs_result *prepare_result(dgeigsf_data *, eigs_result *);


// Eigenvalues and eigenvectors
void dgeigsf(a_int n,
             deigs_phi *phi,
             void *phi_data,
             bool evs,
             const char *which,
             a_int k,
             double tol,
             a_int maxiter,
             eigs_result *result) {

    // Initialize data
    dgeigsf_data *data = dgeigsf_init(n,
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
    dgeigsf_data_destroy(data);
}

// Initialize eigenproblem
static dgeigsf_data *dgeigsf_init(a_int n,
                                  deigs_phi *phi,
                                  void *phi_data,
                                  a_int k,
                                  const char *which,
                                  bool evs,
                                  double tol,
                                  a_int maxiter) {

    // Allocate memory for data
    dgeigsf_data *data = (dgeigsf_data *)malloc(sizeof(dgeigsf_data));

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
    data->resid = (double *)malloc(n*sizeof(double));
    if ((data->ncv = 2*k+1) < 20) data->ncv = 20;
    if (data->ncv > n) data->ncv = n;
    data->v = (double *)calloc(n*data->ncv, sizeof(double));
    data->ldv = n;
    data->iparam = (a_int *)calloc(11, sizeof(a_int));
    data->iparam[0] = 1;
    data->iparam[2] = maxiter;
    data->iparam[3] = 1;
    data->iparam[6] = 1;
    data->ipntr = (a_int *)calloc(14, sizeof(a_int));
    data->workd = (double *)calloc(3*n, sizeof(double));
    data->lworkl = 3*data->ncv*(data->ncv+2);
    data->workl = (double *)calloc(data->lworkl, sizeof(double));
    data->info = 0;
    data->ldz = n;
    data->workev = (double *)calloc(3*data->ncv, sizeof(double));

    // Results
    data->dr = (double *)calloc(data->nev+1, sizeof(double));
    data->di = (double *)calloc(data->nev+1, sizeof(double));
    data->z = (double *)calloc(n*(data->nev+1), sizeof(double));

    return data;
}

// Free for zeigsf_data type
static void dgeigsf_data_destroy(dgeigsf_data *data) {
    free(data->resid); data->resid = NULL;
    free(data->v); data->v = NULL;
    free(data->iparam); data->iparam = NULL;
    free(data->ipntr); data->ipntr = NULL;
    free(data->workd); data->workd = NULL;
    free(data->workl); data->workl = NULL;
    free(data->workev); data->workev = NULL;
    free(data->dr); data->dr = NULL;
    free(data->di); data->di = NULL;
    free(data->z); data->z = NULL;
    free(data);
}

// Do Arnoldi iterations
static void arnoldi_iterations(dgeigsf_data *data) {

    // Arnoldi iterations
    do {
        iterate(data);
    } while ((data->ido == 1) || (data->ido == -1));

    // Check for errors
    if (data->ido != 99) {
        printf("%s\n", "DEIGSF: ARNOLDI PROCESS DID NOT CONVERGE");
        exit(1);
    }
}

// Do a single Arnoldi iteration
static void iterate(dgeigsf_data *data) {

    // Call ZNAUPD
    dnaupd_c(&data->ido,
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
             &data->info);

    // Check for errors
    int nerror = 0;
    if ((data->ido != 1) && (data->ido != -1) && (data->ido != 99)) {
        printf("DEIGSF: ERROR DURING ITERATION: IDO = %d\n", data->ido);
        nerror++;
    }
    if ((data->info != 0) && (data->info != 1)) {
        printf("DEIGSF: ERROR DURING ITERATION: INFO = %d\n", data->info);
        nerror++;
    }
    if (data->info == 1) {
        printf("%s\n", "DEIGSF: MAXIMAL ALLOWED ITERATIONS REACHED");
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
static void extract(dgeigsf_data *data) {

    // For internal use
    const char *howmny = "A"; if (data->evs) howmny = "P";
    a_int *select = (a_int *)calloc(data->ncv, sizeof(a_int));
    double sigmar, sigmai; sigmar = sigmai = 0.; // Not referenced

    // Call ZNEUPD // HERE THINGS GET FUCKED UP            aaaaaaaaaaaaaaaaaaaaaaaaaaaaa
    printf("%d\n", 1000000);
    dneupd_c(data->evs,
             howmny,
             select,
             data->dr,
             data->di,
             data->z,
             data->ldz,
             sigmar,
             sigmai,
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
             &data->info);

    // Clean up
    free(select);

    // Check for errors
    if (data->info) {
        printf("DEIGSF: COULD NOT EXTRACT RESULTS: INFO = %d\n", data->info);
        exit(1);
    }
}

// Load data into result and reorder it to row major
static eigs_result *prepare_result(dgeigsf_data *data, eigs_result *result) {
    a_int n, k, i, j, count;
    n = data->n; k = data->nev; count = 0;
    result->n = n; result->k = k;
    for (j=0; j<k; j++) result->eigvals[j] = CMPLX(data->dr[j], data->di[j]);
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
