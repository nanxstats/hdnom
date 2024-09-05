/*
 *  survAUC_SongZhou.c
 *  SurvM
 *
 *  Created by Sergej Potapov on 14.10.10.
 *  Copyright 2010 __IMBE__. All rights reserved.
 *
 */

#define STRICT_R_HEADERS

#include <Rinternals.h>
#include <Rdefines.h>
#include <R.h>
#include <Rmath.h>

#include "utils.h"

SEXP auc_SZ(SEXP THRESH, SEXP T, SEXP STIME, SEXP EVENT, SEXP N_TIME,
			SEXP STIME_NEW, SEXP EVENT_NEW, SEXP N_TIME_NEW,
			SEXP LP, SEXP N_LP, SEXP LPNEW, SEXP N_LPNEW, SEXP TYPE_SENS)
{
	int nrx, ncx, i, j, k;
	SEXP S1a, xdims, spec;

	double *lp_new;
	lp_new = R_Calloc(INTEGER(N_LPNEW)[0], double);
	for (i = 0; i < INTEGER(N_LPNEW)[0]; i++)
	{
		lp_new[i] = REAL(LPNEW)[i];
	}

	PROTECT(S1a = survfit_cox(LP, STIME, EVENT, N_TIME, N_LP, LPNEW, N_LPNEW));
	xdims = getAttrib(VECTOR_ELT(S1a, 0), R_DimSymbol);
	nrx = INTEGER(xdims)[0];
	ncx = INTEGER(xdims)[1];

	int N_times = LENGTH(T);
	double *surv_new;
	surv_new = R_Calloc(N_times * ncx, double);

	step_eval3(surv_new, REAL(T), REAL(VECTOR_ELT(S1a, 0)), REAL(VECTOR_ELT(S1a, 1)), N_times, ncx, nrx);
	UNPROTECT(1);

	int N_TH = LENGTH(THRESH);
	int N_lpnew = INTEGER(N_LPNEW)[0];
	PROTECT(spec = allocMatrix(REALSXP, N_times, N_TH + 1));
	for (i = 0; i < N_times; i++)
	{
		REAL(spec)
		[i] = 0.0;
	}
	/* Calculation of specificity */
	double tmp_sens_z, tmp_sens_n;
	for (i = 1; i < N_TH + 1; i++)
	{
		for (j = 0; j < N_times; j++)
		{
			tmp_sens_z = 0.;
			tmp_sens_n = 0.;
			for (k = 0; k < N_lpnew; k++)
			{
				if (lp_new[k] > REAL(THRESH)[i - 1])
				{
					tmp_sens_z += surv_new[k * N_times + j];
				}
				tmp_sens_n += surv_new[k * N_times + j];
			}
			REAL(spec)
			[j + N_times * i] = 1.0 - (tmp_sens_z / tmp_sens_n);
		}
	}

	/* Calculation of sensetivity */
	SEXP sens;
	PROTECT(sens = allocMatrix(REALSXP, N_times, N_TH + 1));
	/* last value of sensetivity */
	for (i = 0; i < N_times; i++)
	{
		REAL(sens)
		[i] = 1.0;
	}
	/* type_sens = 0: incident
	   type_sens = 1: cumulative */
	if (!LOGICAL(TYPE_SENS)[0])
	{
		for (i = 1; i < N_TH + 1; i++)
		{
			for (j = 0; j < N_times; j++)
			{
				tmp_sens_z = 0.;
				tmp_sens_n = 0.;
				for (k = 0; k < N_lpnew; k++)
				{
					if (lp_new[k] > REAL(THRESH)[i - 1])
					{
						tmp_sens_z += exp(lp_new[k]) * surv_new[k * N_times + j];
					}
					tmp_sens_n += exp(lp_new[k]) * surv_new[k * N_times + j];
				}
				if (R_finite(tmp_sens_n) && tmp_sens_n > FLT_EPSILON)
				{
					REAL(sens)
					[j + N_times * i] = tmp_sens_z / tmp_sens_n;
				}
				else
				{
					REAL(sens)
					[j + N_times * i] = 0.0;
				}
			}
		}
	}
	else
	{
		double tmp_sens_z, tmp_sens_n;
		for (i = 1; i < N_TH + 1; i++)
		{
			for (j = 0; j < N_times; j++)
			{
				tmp_sens_z = 0.;
				tmp_sens_n = 0.;
				for (k = 0; k < N_lpnew; k++)
				{
					if (lp_new[k] > REAL(THRESH)[i - 1])
					{
						tmp_sens_z += 1.0 - surv_new[k * N_times + j];
					}
					tmp_sens_n += 1.0 - surv_new[k * N_times + j];
				}
				if (R_finite(tmp_sens_n) && tmp_sens_n > FLT_EPSILON)
				{
					REAL(sens)
					[j + N_times * i] = tmp_sens_z / tmp_sens_n;
				}
				else
				{
					REAL(sens)
					[j + N_times * i] = 0.0;
				}
			}
		}
	}

	R_Free(lp_new);
	R_Free(surv_new);
	/* Calculation of AUC */
	SEXP AUC;
	PROTECT(AUC = allocVector(REALSXP, N_times));
	for (i = 0; i < N_times; i++)
	{
		REAL(AUC)
		[i] = 0.;
		for (j = 0; j < N_TH; j++)
		{
			REAL(AUC)
			[i] += ((REAL(sens)[i + N_times * j] + REAL(sens)[i + N_times * (1 + j)]) / 2.0) * fabs((1.0 - REAL(spec)[i + N_times * j]) - (1.0 - REAL(spec)[i + N_times * (1 + j)]));
		}
	}
	SEXP IAUC;
	PROTECT(IAUC = allocVector(REALSXP, 1));
	REAL(IAUC)
	[0] = 0.;

	/* Calculation of iAUC */

	double *f, *S, *S_new;
	int n_new_data = INTEGER(N_TIME_NEW)[0];

	f = R_Calloc(N_times, double);
	S_new = R_Calloc(n_new_data, double);
	S = R_Calloc(N_times, double);
	km_Daim(S_new, REAL(STIME_NEW), REAL(EVENT_NEW), INTEGER(N_TIME_NEW));
	step_eval2(S, REAL(T), S_new, REAL(STIME_NEW), N_times, n_new_data);

	f[0] = 1.0 - S[0];
	for (i = 1; i < N_times; i++)
	{
		f[i] = S[i - 1] - S[i];
	}

	/*
	   type_sens = 0: incident
	   type_sens = 1: cumulative
	 */
	if (!LOGICAL(TYPE_SENS)[0])
	{
		/* incident case*/
		double wT = 0.0;
		for (i = 0; i < N_times; i++)
		{
			wT += f[i] * S[i];
		}
		for (i = 0; i < N_times; i++)
		{
			if (wT != 0.0)
			{
				if (f[i] > FLT_EPSILON)
					REAL(IAUC)
				[0] += REAL(AUC)[i] * (f[i] * S[i]) / wT;
			}
		}
	}
	else
	{
		/* cumulative case*/
		double wT = 0.0;
		for (i = 0; i < N_times; i++)
		{
			if (f[i] > FLT_EPSILON)
				wT += f[i];
		}
		for (i = 0; i < N_times; i++)
		{
			if (wT != 0.0)
			{
				if (f[i] > FLT_EPSILON)
					REAL(IAUC)
				[0] += REAL(AUC)[i] * (f[i]) / wT;
			}
		}
	}
	R_Free(f);
	R_Free(S);
	R_Free(S_new);

	SEXP result, names_result;
	PROTECT(result = allocVector(VECSXP, 5));
	PROTECT(names_result = allocVector(STRSXP, 5));
	SET_STRING_ELT(names_result, 0, mkChar("auc"));
	SET_STRING_ELT(names_result, 1, mkChar("times"));
	SET_STRING_ELT(names_result, 2, mkChar("sens"));
	SET_STRING_ELT(names_result, 3, mkChar("spec"));
	SET_STRING_ELT(names_result, 4, mkChar("iauc"));
	setAttrib(result, R_NamesSymbol, names_result);
	UNPROTECT(1);
	SET_VECTOR_ELT(result, 0, AUC);
	SET_VECTOR_ELT(result, 1, T);
	SET_VECTOR_ELT(result, 2, sens);
	SET_VECTOR_ELT(result, 3, spec);
	SET_VECTOR_ELT(result, 4, IAUC);

	UNPROTECT(5);
	return (result);
}
