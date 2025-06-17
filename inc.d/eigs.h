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


#ifndef EIGS_H
#define EIGS_H

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>
#include <complex.h>
#include <string.h>

#include "../ARPACK/ICB.D/arpack.h"
#include <lapacke.h>


typedef void zeigs_phi(void *,
                       int32_t,
                       const double complex *,
                       double complex *);
typedef void deigs_phi(void *,
                       int32_t,
                       const double *,
                       double *);

typedef struct _EigsResult {
    int32_t n;
    int32_t k;
    double complex *eigvals;
    double complex *eigvecs;
} eigs_result;


eigs_result *eigs(const char *,
                  zeigs_phi *,
                  deigs_phi *,
                  const double complex *,
                  const double *,
                  void *,
                  int32_t,
                  int32_t,
                  const char *,
                  int32_t,
                  double,
                  bool);

void eigs_result_free(eigs_result *);


/* --- Solvers for internal usage ------------------------------------------- */
void zgeigsf(a_int,
             zeigs_phi *,
             void *,
             bool,
             const char *,
             a_int,
             double,
             a_int,
             eigs_result *);
void dgeigsf(a_int,
             deigs_phi *,
             void *,
             bool,
             const char *,
             a_int,
             double,
             a_int,
             eigs_result *);
void zgeigsa(uint32_t,
             const double complex *,
             bool,
             eigs_result *);
void dgeigsa(uint32_t,
             const double *,
             bool,
             eigs_result *);
void zheigsa(uint32_t,
             const double complex *,
             bool,
             eigs_result *);
void dseigsa(uint32_t,
             const double *,
             bool,
             eigs_result *);
/* -------------------------------------------------------------------------- */

#endif
