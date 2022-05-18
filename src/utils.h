/*
 *  utils.h
 *  survM
 *
 *  Created by Sergej Potapov on 07.06.10.
 *  Copyright 2010 __IMBE__. All rights reserved.
 *
 */

#define USE_FC_LEN_T
#include <Rconfig.h>
#include <R_ext/BLAS.h>
#ifndef FCONE
#define FCONE
#endif

#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <R_ext/Applic.h> /* for dgemm */
#include <R_ext/Lapack.h>

/* teked from pnmath Package*/
/*							*/

/* Private header file for use during compilation of Mathlib */
#ifndef MATHLIB_PRIVATE_H
#define MATHLIB_PRIVATE_H

#ifndef MATHLIB_STANDALONE
/* Mathlib in R */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#if defined(HAVE_GLIBC2) && !defined(_BSD_SOURCE)
/* ensure isnan is declared */
#define _BSD_SOURCE 1
#endif
#if defined(HAVE_GLIBC2) && !defined(_ISOC99_SOURCE)
/* ensure expm1 and log1p are declared */
#define _ISOC99_SOURCE 1
#endif
#endif

/* need to add LDOUBLE to Rconfig.h and include Rconfig.h in pnmath.h */
#ifdef DODO
#ifdef HAVE_LONG_DOUBLE
#define LDOUBLE long double
#else
#define LDOUBLE double
#endif
#else
#define LDOUBLE long double
#endif
#endif

/* standalone */
double dmax(double *X, int n);
double dmin(double *X, int n);
double d_mean(double *X, int n);
void step_eval_R(double *s_new, double *t_new, double *s, double *t, int *n_new, int *n);
void step_eval2(double *s_new, double *t_new, double *s, double *t, int n_new, int n);
void step_eval2_left(double *s_new, double *t_new, double *s, double *t, int n_new, int n);
void step_eval3(double *s_new, double *t_new, double *s, double *t, int n_new, int n_s, int n_t);
void rsort_with_x(double *x, double *indx, int n);
void rsort_index(double *x, int *indx, int n);
void rsort_xyz(double *x, double *y, double *indx, int n);
void rsort_xyzv(double *x, double *y, double *z, double *indx, int n);
void km_weight(double *surv, double *time, double *status, double *wt, double *entry, int *n_time);
void km_Daim(double *surv, double *time, double *status, int *n_time);
void cum_sum(double *x, int size);
void my_rev_d(double *x, int *n_x);
void survM_tcrossprod(double *x, int nrx, int ncx, double *y, int nry, int ncy, double *z);
void My_matprod(double *x, int nrx, int ncx, double *y, int nry, int ncy, double *z);
SEXP survfit_cox(SEXP LP, SEXP TIME, SEXP EVENT, SEXP N_TIME, SEXP N_LP, SEXP LPNEW, SEXP N_LPNEW);
