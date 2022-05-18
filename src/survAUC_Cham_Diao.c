/*
 *  survAUC_Cham_Diao.c
 *  Daim
 *
 *  Created by Sergej Potapov on 10.10.10.
 *  Copyright 2010 __IMBE__. All rights reserved.
 *
 */

#include <Rinternals.h>
#include <Rdefines.h>
#include <R.h>
#include <Rmath.h>

#include "utils.h"

SEXP Cham_Diao(SEXP LP, SEXP TH_TIME, SEXP TIME, SEXP EVENT, SEXP N_TIME,
			   SEXP TIME_NEW, SEXP EVENT_NEW, SEXP N_TIME_NEW,
			   SEXP N_LP, SEXP LPNEW, SEXP N_LPNEW)
{

	SEXP xdims, surv, time, S1a;
	int i, j, k, N_th_times, nrx, ncx;

	PROTECT(S1a = survfit_cox(LP, TIME, EVENT, N_TIME, N_LP, LPNEW, N_LPNEW));
	xdims = getAttrib(VECTOR_ELT(S1a, 0), R_DimSymbol);
	nrx = INTEGER(xdims)[0];
	ncx = INTEGER(xdims)[1];
	PROTECT(time = VECTOR_ELT(S1a, 1));
	PROTECT(surv = VECTOR_ELT(S1a, 0));

	N_th_times = LENGTH(TH_TIME);
	double *surv_new;
	surv_new = Calloc(N_th_times * ncx, double);

	step_eval3(surv_new, REAL(TH_TIME), REAL(VECTOR_ELT(S1a, 0)), REAL(VECTOR_ELT(S1a, 1)), N_th_times, ncx, nrx);
	UNPROTECT(1);

	double *factor1;
	factor1 = Calloc(N_th_times, double);

	for (i = 0; i < N_th_times; i++)
	{
		double sumf = 0.0;
		for (j = 0; j < ncx; j++)
		{
			sumf += surv_new[i + j * N_th_times];
		}
		sumf /= (double)ncx;
		factor1[i] = 1.0 / ((1.0 - sumf) * sumf);
	}

	double *EW;
	EW = Calloc(N_th_times, double);

	int n_lpnew = INTEGER(N_LPNEW)[0];
	for (i = 0; i < n_lpnew; i++)
	{
		for (j = 0; j < n_lpnew; j++)
		{
			if (REAL(LPNEW)[j] < REAL(LPNEW)[i])
			{
				for (k = 0; k < N_th_times; k++)
				{
					EW[k] += (1 - surv_new[k + i * N_th_times]) * surv_new[k + j * N_th_times];
				}
			}
		}
	}
	SEXP AUC, IAUC, result, names_result;
	PROTECT(AUC = allocVector(REALSXP, N_th_times));
	for (i = 0; i < N_th_times; i++)
	{
		if (R_finite(factor1[i]))
		{
			REAL(AUC)
			[i] = (factor1[i] * EW[i]) / pow(n_lpnew, 2);
		}
		else
		{
			REAL(AUC)
			[i] = 0.0;
		}
	}
	Free(EW);
	Free(factor1);
	Free(surv_new);
	PROTECT(IAUC = allocVector(REALSXP, 1));
	if (TIME_NEW == R_NilValue)
	{
		REAL(IAUC)
		[0] = 0.0;
	}
	else
	{
		/* Calculation of iAUC */
		int n_new_data = INTEGER(N_TIME_NEW)[0];
		double *f, *S, *S_new;
		f = Calloc(N_th_times, double);
		S_new = Calloc(n_new_data, double);
		S = Calloc(N_th_times, double);
		km_Daim(S_new, REAL(TIME_NEW), REAL(EVENT_NEW), INTEGER(N_TIME_NEW));
		step_eval2(S, REAL(TH_TIME), S_new, REAL(TIME_NEW), N_th_times, n_new_data);

		f[0] = 1.0 - S[0];
		for (i = 1; i < N_th_times; i++)
		{
			f[i] = S[i - 1] - S[i];
		}
		double wT = 0.0;
		for (i = 0; i < N_th_times; i++)
		{
			if (f[i] > FLT_EPSILON)
				wT += f[i];
		}
		double i_auc = 0.0;
		for (i = 0; i < N_th_times; i++)
		{
			/* cumulative case*/
			if (wT != 0.0)
			{
				if (f[i] > FLT_EPSILON && R_finite(REAL(AUC)[i]))
					i_auc += REAL(AUC)[i] * (f[i]) / wT;
			}
		}
		Free(f);
		Free(S);
		Free(S_new);
		REAL(IAUC)
		[0] = i_auc;
	}

	PROTECT(result = allocVector(VECSXP, 3));
	PROTECT(names_result = allocVector(STRSXP, 3));
	SET_STRING_ELT(names_result, 0, mkChar("auc"));
	SET_STRING_ELT(names_result, 1, mkChar("times"));
	SET_STRING_ELT(names_result, 2, mkChar("iauc"));
	setAttrib(result, R_NamesSymbol, names_result);
	UNPROTECT(1);
	SET_VECTOR_ELT(result, 0, AUC);
	SET_VECTOR_ELT(result, 1, TH_TIME);
	SET_VECTOR_ELT(result, 2, IAUC);

	UNPROTECT(5);
	return (result);
}
