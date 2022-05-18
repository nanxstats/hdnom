/*
 *  survAUC_UNO.h
 *  Daim
 *
 *  Created by Sergej Potapov on 01.06.10.
 *  Copyright 2010 __IMBE__. All rights reserved.
 *
 */

void auc_uno(double *auc, double *i_auc, double *sens, double *spec, double *surv_time,
			 double *status, double *thres, double *t, double *marker, double *new_surv_t,
			 double *new_event, int *n_th, int *n_t, int *n_new_data, int *n_surv);
