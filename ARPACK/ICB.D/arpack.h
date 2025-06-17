#ifndef ARPACK_H
#define ARPACK_H

#include <stdint.h>
#include <complex.h>


#define a_int int32_t
#define a_dcomplex double _Complex

#ifndef CMPLX
#define CMPLX(r,i) ((double _Complex)((double)(r) + _Complex_I * (double)(i)))
#endif


// Double complex routines for general and hermitian endomorphisms
void znaupd_c(a_int*            ido      ,
              char const*       bmat     ,
              a_int             n        ,
              char const*       which    ,
              a_int             nev      ,
              double            tol      ,
              a_dcomplex*       resid    ,
              a_int             ncv      ,
              a_dcomplex*       v        ,
              a_int             ldv      ,
              a_int*            iparam   ,
              a_int*            ipntr    ,
              a_dcomplex*       workd    ,
              a_dcomplex*       workl    ,
              a_int             lworkl   ,
              double*           rwork    ,
              a_int*            info      );
void zneupd_c(a_int             rvec     ,
              char const*       howmny   ,
              a_int const*      select   ,
              a_dcomplex*       d        ,
              a_dcomplex*       z        ,
              a_int             ldz      ,
              a_dcomplex        sigma    ,
              a_dcomplex*       workev   ,
              char const*       bmat     ,
              a_int             n        ,
              char const*       which    ,
              a_int             nev      ,
              double            tol      ,
              a_dcomplex*       resid    ,
              a_int             ncv      ,
              a_dcomplex*       v        ,
              a_int             ldv      ,
              a_int*            iparam   ,
              a_int*            ipntr    ,
              a_dcomplex*       workd    ,
              a_dcomplex*       workl    ,
              a_int             lworkl   ,
              double*           rwork    ,
              a_int*            info      );

// Double routines for general endomorphisms
void dnaupd_c(a_int*            ido      ,
              char const*       bmat     ,
              a_int             n        ,
              char const*       which    ,
              a_int             nev      ,
              double            tol      ,
              double*           resid    ,
              a_int             ncv      ,
              double*           v        ,
              a_int             ldv      ,
              a_int*            iparam   ,
              a_int*            ipntr    ,
              double*           workd    ,
              double*           workl    ,
              a_int             lworkl   ,
              a_int*            info      );
void dneupd_c(a_int             rvec     ,
              char const*       howmny   ,
              a_int const*      select   ,
              double*           dr       ,
              double*           di       ,
              double*           z        ,
              a_int             ldz      ,
              double            sigmar   ,
              double            sigmai   ,
              double*           workev   ,
              char const*       bmat     ,
              a_int             n        ,
              char const*       which    ,
              a_int             nev      ,
              double            tol      ,
              double*           resid    ,
              a_int             ncv      ,
              double*           v        ,
              a_int             ldv      ,
              a_int*            iparam   ,
              a_int*            ipntr    ,
              double*           workd    ,
              double*           workl    ,
              a_int             lworkl   ,
              a_int*            info      );

#endif
