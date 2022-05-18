/*
 *  utils.c
 *  survM
 *
 *  Created by Sergej Potapov on 07.06.10.
 *  Copyright 2010 __IMBE__. All rights reserved.
 *
 */

#include "utils.h"

double dmax(double *X, int n)
{
	int i = 0;
	double max = 0.0;
	for (i = 0; i < n; i++)
	{
		if (X[i] > max)
		{
			max = X[i];
		}
	}
	return (max);
}

double dmin(double *X, int n)
{
	int i = 0;
	double min = 0.0;
	for (i = 0; i < n; i++)
	{
		if (X[i] < min)
		{
			min = X[i];
		}
	}
	return (min);
}

double d_mean(double *X, int n)
{
	int i = 0;
	double max = 0.0;
	for (i = 0; i < n; i++)
	{
		max += X[i];
	}
	max /= (double)n;
	return (max);
}

/*	help function for evaluation of the survival step-function
 *\param S.new vector
 *\param t.new vector of the new time points, which will be evaluated
 *\param s - vector of the suvival
 *\param t - vector of the time points
 *\param n.new - length of the vectors S.new and t.new
 *\param n - length of the vectors S and t.
 */

void step_eval_R(double *s_new, double *t_new, double *s, double *t, int *n_new, int *n)
{
	int i, j, optim;

	for (j = 0; j < *n_new; j++)
	{
		optim = 1;
		for (i = *n - 1; i >= 0; i--)
		{
			if (optim && t_new[j] >= t[i])
			{
				s_new[j] = s[i];
				optim = 0;
			}
		}
		if (optim)
		{
			s_new[j] = 1;
		}
	}
}

/*	help function for evaluation of the survival step-function
 *\param S.new vector
 *\param t.new vector of the new time points, which will be evaluated
 *\param s - vector of the suvival
 *\param t - vector of the time points
 *\param n.new - length of the vectors S.new and t.new
 *\param n - length of the vectors S and t.
 */

void step_eval2(double *s_new, double *t_new, double *s, double *t, int n_new, int n)
{
	int i, j, optim;

	for (j = 0; j < n_new; j++)
	{
		optim = 1;
		for (i = n - 1; i >= 0; i--)
		{
			if (optim && t_new[j] >= t[i])
			{
				s_new[j] = s[i];
				optim = 0;
			}
		}
		if (optim)
		{
			s_new[j] = 1;
		}
	}
}

/* left-sided limit */

void step_eval2_left(double *s_new, double *t_new, double *s, double *t, int n_new, int n)
{
	int i, j, optim;

	for (j = 0; j < n_new; j++)
	{
		optim = 1;
		for (i = n - 1; i >= 0; i--)
		{
			if (optim && (t_new[j] - FLT_EPSILON) >= t[i])
			{
				s_new[j] = s[i];
				optim = 0;
			}
		}
		if (optim)
		{
			s_new[j] = 1;
		}
	}
}

/*	help function for evaluation of the survival step-function
 *\param S.new vector
 *\param t.new vector of the new time points, which will be evaluated
 *\param s - matrix of the suvival
 *\param t - vector of the time points
 *\param n.new - length of the vectors S.new and t.new
 *\param n_s - number of cols of the matrix S
 *\param n_t - length of vector t
 */

void step_eval3(double *s_new, double *t_new, double *s, double *t, int n_new, int n_s, int n_t)
{
	int i, j, k, optim;

	for (k = 0; k < n_s; k++)
	{
		for (j = 0; j < n_new; j++)
		{
			optim = 1;
			for (i = n_t - 1; i >= 0; i--)
			{
				if (optim && t_new[j] >= t[i])
				{
					s_new[j + n_new * k] = s[i + n_t * k];
					optim = 0;
				}
			}
			if (optim)
			{
				s_new[j + n_new * k] = 1;
			}
		}
	}
}

static int rcmp_TW(double x, double y, Rboolean nalast)
{
	int nax = ISNAN(x), nay = ISNAN(y);
	if (nax && nay)
		return 0;
	if (nax)
		return nalast ? 1 : -1;
	if (nay)
		return nalast ? -1 : 1;
	if (x < y)
		return -1;
	if (x > y)
		return 1;
	return 0;
}

void rsort_with_x(double *x, double *indx, int n)
{
	double v, iv;
	int i, j, h;

	for (h = 1; h <= n / 9; h = 3 * h + 1)
		;
	for (; h > 0; h /= 3)
		for (i = h; i < n; i++)
		{
			v = x[i];
			iv = indx[i];
			j = i;
			while (j >= h && rcmp_TW(x[j - h], v, TRUE) > 0)
			{
				x[j] = x[j - h];
				indx[j] = indx[j - h];
				j -= h;
			}
			x[j] = v;
			indx[j] = iv;
		}
}

void rsort_index(double *x, int *indx, int n)
{
	double v;
	int i, j, h, iv;

	for (h = 1; h <= n / 9; h = 3 * h + 1)
		;
	for (; h > 0; h /= 3)
		for (i = h; i < n; i++)
		{
			v = x[i];
			iv = indx[i];
			j = i;
			while (j >= h && rcmp_TW(x[j - h], v, TRUE) > 0)
			{
				x[j] = x[j - h];
				indx[j] = indx[j - h];
				j -= h;
			}
			x[j] = v;
			indx[j] = iv;
		}
}

void rsort_xyz(double *x, double *y, double *indx, int n)
{
	double v, iv, vi;
	int i, j, h;

	for (h = 1; h <= n / 9; h = 3 * h + 1)
		;
	for (; h > 0; h /= 3)
		for (i = h; i < n; i++)
		{
			v = x[i];
			iv = indx[i];
			vi = y[i];
			j = i;
			while (j >= h && rcmp_TW(x[j - h], v, TRUE) > 0)
			{
				x[j] = x[j - h];
				indx[j] = indx[j - h];
				y[j] = y[j - h];
				j -= h;
			}
			x[j] = v;
			indx[j] = iv;
			y[j] = vi;
		}
}

void rsort_xyzv(double *x, double *y, double *z, double *indx, int n)
{
	double v, iv, vi, vii;
	int i, j, h;

	for (h = 1; h <= n / 9; h = 3 * h + 1)
		;
	for (; h > 0; h /= 3)
		for (i = h; i < n; i++)
		{
			v = x[i];
			iv = indx[i];
			vi = y[i];
			vii = z[i];
			j = i;
			while (j >= h && rcmp_TW(x[j - h], v, TRUE) > 0)
			{
				x[j] = x[j - h];
				indx[j] = indx[j - h];
				y[j] = y[j - h];
				z[j] = z[j - h];
				j -= h;
			}
			x[j] = v;
			indx[j] = iv;
			y[j] = vi;
			z[j] = vii;
		}
}

/* weighted KM - from risksetROC */

void km_weight(double *surv, double *time, double *status, double *wt, double *entry, int *n_time)
{
	int i, j, dead, at_risk;
	rsort_with_x(time, status, *n_time);

	double current = 1.0;
	for (i = 0; i < *n_time; i++)
	{
		dead = 0;
		at_risk = 0;
		for (j = 0; j < *n_time; j++)
		{
			at_risk += (entry[i] <= time[j]) && (time[i] <= time[j]) && wt[i];
			dead += (entry[i] <= time[j]) && (time[i] == time[j]) && (status[i]) && wt[i];
		}
		current = current * (1.0 - (double)dead / (double)at_risk);
		surv[i] = current;
	}
}

/* Kaplan-Meier estimation  */

void km_Daim(double *surv, double *time, double *status, int *n_time)
{
	int i, j, dead, at_risk;
	rsort_with_x(time, status, *n_time);

	double current = 1.0;
	for (i = 0; i < *n_time; i++)
	{
		dead = 0;
		at_risk = 0;
		for (j = 0; j < *n_time; j++)
		{
			at_risk += (time[i] <= time[j]);
			dead += (time[i] == time[j]) && (status[i]);
		}
		current = current * (1.0 - (double)dead / (double)at_risk);
		surv[i] = current;
	}
}

void cum_sum(double *x, int size)
{
	LDOUBLE sum = 0.;
	int i;
	for (i = 0; i < size; i++)
	{
		sum += x[i];
		x[i] = sum;
	}
}

void my_rev_d(double *x, int *n_x)
{
	double swap = 0.0;
	int i, j = *n_x - 1;
	for (i = 0; i < j; i++, j--)
	{
		swap = x[i];
		x[i] = x[j];
		x[j] = swap;
	}
}

void My_matprod(double *x, int nrx, int ncx, double *y, int nry, int ncy, double *z)
{
	char *transa = "N", *transb = "N";
	int i, j, k;
	double one = 1.0, zero = 0.0;
	LDOUBLE sum;
	Rboolean have_na = FALSE;

	if (nrx > 0 && ncx > 0 && nry > 0 && ncy > 0)
	{
		/* Don't trust the BLAS to handle NA/NaNs correctly: PR#4582
		 * The test is only O(n) here
		 */
		for (i = 0; i < nrx * ncx; i++)
			if (ISNAN(x[i]))
			{
				have_na = TRUE;
				break;
			}
		if (!have_na)
			for (i = 0; i < nry * ncy; i++)
				if (ISNAN(y[i]))
				{
					have_na = TRUE;
					break;
				}
		if (have_na)
		{
			for (i = 0; i < nrx; i++)
				for (k = 0; k < ncy; k++)
				{
					sum = 0.0;
					for (j = 0; j < ncx; j++)
						sum += x[i + j * nrx] * y[j + k * nry];
					z[i + k * nrx] = sum;
				}
		}
		else
			F77_CALL(dgemm)
			(transa, transb, &nrx, &ncy, &ncx, &one,
			 x, &nrx, y, &nry, &zero, z, &nrx FCONE FCONE);
	}
	else /* zero-extent operations should return zeroes */
		for (i = 0; i < nrx * ncy; i++)
			z[i] = 0;
}

void survM_tcrossprod(double *x, int nrx, int ncx,
					  double *y, int nry, int ncy, double *z)
{
	char *transa = "N", *transb = "T";
	double one = 1.0, zero = 0.0;
	if (nrx > 0 && ncx > 0 && nry > 0 && ncy > 0)
	{
		F77_CALL(dgemm)
		(transa, transb, &nrx, &nry, &ncx, &one,
		 x, &nrx, y, &nry, &zero, z, &nrx FCONE FCONE);
	}
	else
	{ /* zero-extent operations should return zeroes */
		int i;
		for (i = 0; i < nrx * nry; i++)
			z[i] = 0;
	}
}

SEXP survfit_cox(SEXP LP, SEXP TIME, SEXP EVENT, SEXP N_TIME, SEXP N_LP, SEXP LPNEW, SEXP N_LPNEW)
{
	int *n_time = INTEGER(N_TIME);
	int *n_lp = INTEGER(N_LP);
	int *n_lpnew = INTEGER(N_LPNEW);
	int i, j, k, time0;

	double *lp, *time, *event, *lpnew;
	time = Calloc(*n_time, double);
	event = Calloc(*n_time, double);
	lp = Calloc(*n_lp, double);
	lpnew = Calloc(*n_lpnew, double);

	for (i = 0; i < *n_lp; i++)
	{
		lp[i] = REAL(LP)[i];
	}
	for (i = 0; i < *n_time; i++)
	{
		time[i] = REAL(TIME)[i];
		event[i] = REAL(EVENT)[i];
	}
	for (i = 0; i < *n_lpnew; i++)
	{
		lpnew[i] = REAL(LPNEW)[i];
	}

	rsort_xyz(time, event, lp, *n_time);

	double *n_event, *R;
	n_event = Calloc(*n_time, double);
	R = Calloc(*n_time, double);

	double lp_mean = 0.0;
	lp_mean = d_mean(lp, *n_lp);

	for (i = 0, j = *n_lp - 1; i < *n_lp; i++, j--)
	{
		lp[i] -= lp_mean;
		R[j] = exp(lp[i]);
	}
	cum_sum(R, *n_lp);
	my_rev_d(R, n_lp);

	int *diff_time;
	diff_time = Calloc(*n_time, int);
	diff_time[0] = 1;

	time0 = time[0];
	n_event[0] = event[0];
	for (i = 1, j = 1; i < *n_time; i++)
	{
		diff_time[i] = 0;
		n_event[j] += event[i];
		if (fabs(time0 - time[i]) > DBL_EPSILON)
		{
			diff_time[i] = 1;
			time0 = time[i];
			j++;
		}
	}
	for (i = 0, k = 0; i < *n_time; i++)
	{
		if (diff_time[i])
		{
			R[k] = n_event[k] / R[i];
			time[k] = time[i];
			k++;
		}
	}
	Free(diff_time);
	cum_sum(R, k);
	for (i = 0, j = 0; i < k; i++)
	{
		if (n_event[i] > 0.)
		{
			R[j] = -1. * R[i];
			time[j] = time[i];
			j++;
		}
	}

	for (i = 0; i < *n_lpnew; i++)
	{
		lpnew[i] = exp(lpnew[i] - lp_mean);
	}
	SEXP ans, utimes, nevent;
	PROTECT(utimes = allocVector(REALSXP, j));
	PROTECT(nevent = allocVector(REALSXP, j));
	PROTECT(ans = allocMatrix(REALSXP, j, *n_lpnew));

	survM_tcrossprod(R, j, 1, lpnew, *n_lpnew, 1, REAL(ans));
	for (i = 0; i < j * (*n_lpnew); i++)
	{
		REAL(ans)
		[i] = exp(REAL(ans)[i]);
	}
	for (i = 0; i < j; i++)
	{
		REAL(utimes)
		[i] = time[i];
		REAL(nevent)
		[i] = n_event[i];
	}
	SEXP result_out = PROTECT(allocVector(VECSXP, 3));

	Free(n_event);
	Free(R);
	Free(lp);
	Free(time);
	Free(event);
	SET_VECTOR_ELT(result_out, 0, ans);
	SET_VECTOR_ELT(result_out, 1, utimes);
	SET_VECTOR_ELT(result_out, 2, nevent);
	UNPROTECT(4);
	return (result_out);
}
