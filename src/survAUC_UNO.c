/*
 *  survAUC_UNO.c
 *  Daim
 *
 *  Created by Sergej Potapov on 01.06.10.
 *  Copyright 2010 __IMBE__. All rights reserved.
 *
 */

#include <Rinternals.h>
#include <Rdefines.h>
#include <R.h>
#include <Rmath.h>

#include "utils.h"

/* Uno AUC
 *\param auc double, vector, which length is length(t)
 *\param sens double, vector, which length is length(thres)*(length(t)+1)
 *\param spec double, vector, which length is length(thres)*(length(t)+1)
 *\param surv double, surival vector
 *\param surv_time double, surival vector
 *\param thres double, vector
 *\param t vector of the time points
 *\param marker - vector, linear predictor
 *\param new_data - vector of the new data
 *\param n_th - length of the threshold vector
 *\param n_t - length of the t vector
 *\param n_new_data - length of the new_data vector
 *\param Con_Inc integer: if Con_Inc=1 wenn Cond, else Dynamic.
 */

void auc_uno(double *auc, double *i_auc, double *sens, double *spec, double *surv_time,
			 double *status, double *thres, double *t, double *marker, double *new_surv_t,
			 double *new_event, int *n_th, int *n_t, int *n_new_data, int *n_surv)
{

	/* Calculation of sensetivity */
	int k, i, j;
	double Ivec_zse, Ivec_nse;

	rsort_with_x(surv_time, status, *n_surv);

	double *SProb;
	SProb = Calloc(*n_surv, double);
	km_Daim(SProb, surv_time, status, n_surv);

	double *G;
	G = Calloc(*n_new_data, double);
	step_eval2(G, new_surv_t, SProb, surv_time, *n_new_data, *n_surv);

	for (k = 1; k < *n_th + 1; k++)
	{
		for (j = 0; j < *n_t; j++)
		{
			Ivec_zse = 0.0, Ivec_nse = 0.0;
			for (i = 0; i < *n_new_data; i++)
			{
				if (t[j] >= new_surv_t[i])
				{
					if (marker[i] > thres[k - 1])
					{
						Ivec_zse += new_event[i] / G[i];
					}
					Ivec_nse += new_event[i] / G[i];
				}
			}
			if (Ivec_nse > FLT_EPSILON)
			{
				sens[k * (*n_t) + j] = Ivec_zse / Ivec_nse;
			}
			else
			{
				sens[k * (*n_t) + j] = 0.0;
			}
		}
	}
	Free(SProb);
	Free(G);
	/* Calculation of specificity */
	double Ivec_zsp, Ivec_nsp, tmp_Ivec_zsp = 0.0;

	for (k = 1; k < *n_th + 1; k++)
	{
		for (j = 0; j < *n_t; j++)
		{
			Ivec_zsp = 0.0, Ivec_nsp = 0.0;
			for (i = 0; i < *n_new_data; i++)
			{
				tmp_Ivec_zsp = t[j] < new_surv_t[i];
				Ivec_zsp += (marker[i] <= thres[k - 1]) * tmp_Ivec_zsp;
				Ivec_nsp += tmp_Ivec_zsp;
			}
			if (Ivec_nsp > FLT_EPSILON)
			{
				spec[k * (*n_t) + j] = Ivec_zsp / Ivec_nsp;
			}
			else
			{
				spec[k * (*n_t) + j] = 0.0;
			}
		}
	}
	/* Calculation of AUC */
	for (i = 0; i < *n_t; i++)
	{
		for (j = 0; j < *n_th; j++)
		{
			auc[i] += ((sens[i + *n_t * j] + sens[i + *n_t * (1 + j)]) / 2.0) * fabs((1.0 - spec[i + *n_t * j]) - (1.0 - spec[i + *n_t * (1 + j)]));
		}
	}
	/* Calculation of iAUC */
	double *f, *S, *S_new;
	f = Calloc(*n_t, double);
	S_new = Calloc(*n_new_data, double);
	S = Calloc(*n_t, double);
	km_Daim(S_new, new_surv_t, new_event, n_new_data);
	step_eval2(S, t, S_new, new_surv_t, *n_t, *n_new_data);

	f[0] = 1.0 - S[0];
	for (i = 1; i < *n_t; i++)
	{
		f[i] = S[i - 1] - S[i];
	}
	double wT = 0.0;
	for (i = 0; i < *n_t; i++)
	{
		if (f[i] > FLT_EPSILON)
		{
			wT += f[i];
		}
	}
	for (i = 0; i < *n_t; i++)
	{
		if (wT != 0.0)
		{
			/* cumulative case*/
			if (f[i] > FLT_EPSILON)
			{
				*i_auc += auc[i] * (f[i]) / wT;
			}
		}
	}
	Free(f);
	Free(S);
	Free(S_new);
}
