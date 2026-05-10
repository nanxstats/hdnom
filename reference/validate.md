# Validate high-dimensional Cox models with time-dependent AUC

Validate high-dimensional Cox models with time-dependent AUC

## Usage

``` r
validate(
  x,
  time,
  event,
  model.type = c("lasso", "alasso", "flasso", "enet", "aenet", "mcp", "mnet", "scad",
    "snet"),
  alpha,
  lambda,
  pen.factor = NULL,
  gamma,
  lambda1,
  lambda2,
  method = c("bootstrap", "cv", "repeated.cv"),
  boot.times = NULL,
  nfolds = NULL,
  rep.times = NULL,
  tauc.type = c("CD", "SZ", "UNO"),
  tauc.time,
  seed = 1001,
  trace = TRUE
)
```

## Arguments

- x:

  Matrix of training data used for fitting the model; on which to run
  the validation.

- time:

  Survival time. Must be of the same length with the number of rows as
  `x`.

- event:

  Status indicator, normally 0 = alive, 1 = dead. Must be of the same
  length with the number of rows as `x`.

- model.type:

  Model type to validate. Could be one of `"lasso"`, `"alasso"`,
  `"flasso"`, `"enet"`, `"aenet"`, `"mcp"`, `"mnet"`, `"scad"`, or
  `"snet"`.

- alpha:

  Value of the elastic-net mixing parameter alpha for `enet`, `aenet`,
  `mnet`, and `snet` models. For `lasso`, `alasso`, `mcp`, and `scad`
  models, please set `alpha = 1`. `alpha=1`: lasso (l1) penalty;
  `alpha=0`: ridge (l2) penalty. Note that for `mnet` and `snet` models,
  `alpha` can be set to very close to 0 but not 0 exactly.

- lambda:

  Value of the penalty parameter lambda to use in the model fits on the
  resampled data. From the fitted Cox model.

- pen.factor:

  Penalty factors to apply to each coefficient. From the fitted
  *adaptive lasso* or *adaptive elastic-net* model.

- gamma:

  Value of the model parameter gamma for MCP/SCAD/Mnet/Snet models.

- lambda1:

  Value of the penalty parameter lambda1 for fused lasso model.

- lambda2:

  Value of the penalty parameter lambda2 for fused lasso model.

- method:

  Validation method. Could be `"bootstrap"`, `"cv"`, or `"repeated.cv"`.

- boot.times:

  Number of repetitions for bootstrap.

- nfolds:

  Number of folds for cross-validation and repeated cross-validation.

- rep.times:

  Number of repeated times for repeated cross-validation.

- tauc.type:

  Type of time-dependent AUC. Including `"CD"` proposed by Chambless and
  Diao (2006)., `"SZ"` proposed by Song and Zhou (2008)., `"UNO"`
  proposed by Uno et al. (2007).

- tauc.time:

  Numeric vector. Time points at which to evaluate the time-dependent
  AUC.

- seed:

  A random seed for resampling.

- trace:

  Logical. Output the validation progress or not. Default is `TRUE`.

## References

Chambless, L. E. and G. Diao (2006). Estimation of time-dependent area
under the ROC curve for long-term risk prediction. *Statistics in
Medicine* 25, 3474–3486.

Song, X. and X.-H. Zhou (2008). A semiparametric approach for the
covariate specific ROC curve with survival outcome. *Statistica Sinica*
18, 947–965.

Uno, H., T. Cai, L. Tian, and L. J. Wei (2007). Evaluating prediction
rules for t-year survivors with censored regression models. *Journal of
the American Statistical Association* 102, 527–537.

## Examples

``` r
data(smart)
x <- as.matrix(smart[, -c(1, 2)])[1:500, ]
time <- smart$TEVENT[1:500]
event <- smart$EVENT[1:500]
y <- survival::Surv(time, event)

fit <- fit_lasso(x, y, nfolds = 5, rule = "lambda.min", seed = 11)
#> Warning: Starting in glmnet 5.1, the default Cox tie-handling method will change from 'breslow' to 'efron' (matching survival::coxph). To silence this message and lock in the v5.0 default, pass cox.ties = 'breslow' explicitly. To preview the v5.1 behavior, pass cox.ties = 'efron'.
#> Warning: Starting in glmnet 5.1, the default Cox tie-handling method will change from 'breslow' to 'efron' (matching survival::coxph). To silence this message and lock in the v5.0 default, pass cox.ties = 'breslow' explicitly. To preview the v5.1 behavior, pass cox.ties = 'efron'.
#> Warning: Starting in glmnet 5.1, the default Cox tie-handling method will change from 'breslow' to 'efron' (matching survival::coxph). To silence this message and lock in the v5.0 default, pass cox.ties = 'breslow' explicitly. To preview the v5.1 behavior, pass cox.ties = 'efron'.
#> Warning: Starting in glmnet 5.1, the default Cox tie-handling method will change from 'breslow' to 'efron' (matching survival::coxph). To silence this message and lock in the v5.0 default, pass cox.ties = 'breslow' explicitly. To preview the v5.1 behavior, pass cox.ties = 'efron'.
#> Warning: Starting in glmnet 5.1, the default Cox tie-handling method will change from 'breslow' to 'efron' (matching survival::coxph). To silence this message and lock in the v5.0 default, pass cox.ties = 'breslow' explicitly. To preview the v5.1 behavior, pass cox.ties = 'efron'.
#> Warning: Starting in glmnet 5.1, the default Cox tie-handling method will change from 'breslow' to 'efron' (matching survival::coxph). To silence this message and lock in the v5.0 default, pass cox.ties = 'breslow' explicitly. To preview the v5.1 behavior, pass cox.ties = 'efron'.
#> Warning: Starting in glmnet 5.1, the default Cox tie-handling method will change from 'breslow' to 'efron' (matching survival::coxph). To silence this message and lock in the v5.0 default, pass cox.ties = 'breslow' explicitly. To preview the v5.1 behavior, pass cox.ties = 'efron'.

# Model validation by bootstrap with time-dependent AUC
# Normally boot.times should be set to 200 or more,
# we set it to 3 here only to save example running time.
val.boot <- validate(
  x, time, event,
  model.type = "lasso",
  alpha = 1, lambda = fit$lambda,
  method = "bootstrap", boot.times = 3,
  tauc.type = "UNO", tauc.time = seq(0.25, 2, 0.25) * 365,
  seed = 1010
)
#> Start bootstrap sample 1 
#> Warning: Starting in glmnet 5.1, the default Cox tie-handling method will change from 'breslow' to 'efron' (matching survival::coxph). To silence this message and lock in the v5.0 default, pass cox.ties = 'breslow' explicitly. To preview the v5.1 behavior, pass cox.ties = 'efron'.
#> Start bootstrap sample 2 
#> Warning: Starting in glmnet 5.1, the default Cox tie-handling method will change from 'breslow' to 'efron' (matching survival::coxph). To silence this message and lock in the v5.0 default, pass cox.ties = 'breslow' explicitly. To preview the v5.1 behavior, pass cox.ties = 'efron'.
#> Start bootstrap sample 3 
#> Warning: Starting in glmnet 5.1, the default Cox tie-handling method will change from 'breslow' to 'efron' (matching survival::coxph). To silence this message and lock in the v5.0 default, pass cox.ties = 'breslow' explicitly. To preview the v5.1 behavior, pass cox.ties = 'efron'.

# Model validation by 5-fold cross-validation with time-dependent AUC
val.cv <- validate(
  x, time, event,
  model.type = "lasso",
  alpha = 1, lambda = fit$lambda,
  method = "cv", nfolds = 5,
  tauc.type = "UNO", tauc.time = seq(0.25, 2, 0.25) * 365,
  seed = 1010
)
#> Start fold 1 
#> Warning: Starting in glmnet 5.1, the default Cox tie-handling method will change from 'breslow' to 'efron' (matching survival::coxph). To silence this message and lock in the v5.0 default, pass cox.ties = 'breslow' explicitly. To preview the v5.1 behavior, pass cox.ties = 'efron'.
#> Start fold 2 
#> Warning: Starting in glmnet 5.1, the default Cox tie-handling method will change from 'breslow' to 'efron' (matching survival::coxph). To silence this message and lock in the v5.0 default, pass cox.ties = 'breslow' explicitly. To preview the v5.1 behavior, pass cox.ties = 'efron'.
#> Start fold 3 
#> Warning: Starting in glmnet 5.1, the default Cox tie-handling method will change from 'breslow' to 'efron' (matching survival::coxph). To silence this message and lock in the v5.0 default, pass cox.ties = 'breslow' explicitly. To preview the v5.1 behavior, pass cox.ties = 'efron'.
#> Start fold 4 
#> Warning: Starting in glmnet 5.1, the default Cox tie-handling method will change from 'breslow' to 'efron' (matching survival::coxph). To silence this message and lock in the v5.0 default, pass cox.ties = 'breslow' explicitly. To preview the v5.1 behavior, pass cox.ties = 'efron'.
#> Start fold 5 
#> Warning: Starting in glmnet 5.1, the default Cox tie-handling method will change from 'breslow' to 'efron' (matching survival::coxph). To silence this message and lock in the v5.0 default, pass cox.ties = 'breslow' explicitly. To preview the v5.1 behavior, pass cox.ties = 'efron'.

# Model validation by repeated cross-validation with time-dependent AUC
val.repcv <- validate(
  x, time, event,
  model.type = "lasso",
  alpha = 1, lambda = fit$lambda,
  method = "repeated.cv", nfolds = 5, rep.times = 3,
  tauc.type = "UNO", tauc.time = seq(0.25, 2, 0.25) * 365,
  seed = 1010
)
#> Start repeat round 1 fold 1 
#> Warning: Starting in glmnet 5.1, the default Cox tie-handling method will change from 'breslow' to 'efron' (matching survival::coxph). To silence this message and lock in the v5.0 default, pass cox.ties = 'breslow' explicitly. To preview the v5.1 behavior, pass cox.ties = 'efron'.
#> Start repeat round 1 fold 2 
#> Warning: Starting in glmnet 5.1, the default Cox tie-handling method will change from 'breslow' to 'efron' (matching survival::coxph). To silence this message and lock in the v5.0 default, pass cox.ties = 'breslow' explicitly. To preview the v5.1 behavior, pass cox.ties = 'efron'.
#> Start repeat round 1 fold 3 
#> Warning: Starting in glmnet 5.1, the default Cox tie-handling method will change from 'breslow' to 'efron' (matching survival::coxph). To silence this message and lock in the v5.0 default, pass cox.ties = 'breslow' explicitly. To preview the v5.1 behavior, pass cox.ties = 'efron'.
#> Start repeat round 1 fold 4 
#> Warning: Starting in glmnet 5.1, the default Cox tie-handling method will change from 'breslow' to 'efron' (matching survival::coxph). To silence this message and lock in the v5.0 default, pass cox.ties = 'breslow' explicitly. To preview the v5.1 behavior, pass cox.ties = 'efron'.
#> Start repeat round 1 fold 5 
#> Warning: Starting in glmnet 5.1, the default Cox tie-handling method will change from 'breslow' to 'efron' (matching survival::coxph). To silence this message and lock in the v5.0 default, pass cox.ties = 'breslow' explicitly. To preview the v5.1 behavior, pass cox.ties = 'efron'.
#> Start repeat round 2 fold 1 
#> Warning: Starting in glmnet 5.1, the default Cox tie-handling method will change from 'breslow' to 'efron' (matching survival::coxph). To silence this message and lock in the v5.0 default, pass cox.ties = 'breslow' explicitly. To preview the v5.1 behavior, pass cox.ties = 'efron'.
#> Start repeat round 2 fold 2 
#> Warning: Starting in glmnet 5.1, the default Cox tie-handling method will change from 'breslow' to 'efron' (matching survival::coxph). To silence this message and lock in the v5.0 default, pass cox.ties = 'breslow' explicitly. To preview the v5.1 behavior, pass cox.ties = 'efron'.
#> Start repeat round 2 fold 3 
#> Warning: Starting in glmnet 5.1, the default Cox tie-handling method will change from 'breslow' to 'efron' (matching survival::coxph). To silence this message and lock in the v5.0 default, pass cox.ties = 'breslow' explicitly. To preview the v5.1 behavior, pass cox.ties = 'efron'.
#> Start repeat round 2 fold 4 
#> Warning: Starting in glmnet 5.1, the default Cox tie-handling method will change from 'breslow' to 'efron' (matching survival::coxph). To silence this message and lock in the v5.0 default, pass cox.ties = 'breslow' explicitly. To preview the v5.1 behavior, pass cox.ties = 'efron'.
#> Start repeat round 2 fold 5 
#> Warning: Starting in glmnet 5.1, the default Cox tie-handling method will change from 'breslow' to 'efron' (matching survival::coxph). To silence this message and lock in the v5.0 default, pass cox.ties = 'breslow' explicitly. To preview the v5.1 behavior, pass cox.ties = 'efron'.
#> Start repeat round 3 fold 1 
#> Warning: Starting in glmnet 5.1, the default Cox tie-handling method will change from 'breslow' to 'efron' (matching survival::coxph). To silence this message and lock in the v5.0 default, pass cox.ties = 'breslow' explicitly. To preview the v5.1 behavior, pass cox.ties = 'efron'.
#> Start repeat round 3 fold 2 
#> Warning: Starting in glmnet 5.1, the default Cox tie-handling method will change from 'breslow' to 'efron' (matching survival::coxph). To silence this message and lock in the v5.0 default, pass cox.ties = 'breslow' explicitly. To preview the v5.1 behavior, pass cox.ties = 'efron'.
#> Start repeat round 3 fold 3 
#> Warning: Starting in glmnet 5.1, the default Cox tie-handling method will change from 'breslow' to 'efron' (matching survival::coxph). To silence this message and lock in the v5.0 default, pass cox.ties = 'breslow' explicitly. To preview the v5.1 behavior, pass cox.ties = 'efron'.
#> Start repeat round 3 fold 4 
#> Warning: Starting in glmnet 5.1, the default Cox tie-handling method will change from 'breslow' to 'efron' (matching survival::coxph). To silence this message and lock in the v5.0 default, pass cox.ties = 'breslow' explicitly. To preview the v5.1 behavior, pass cox.ties = 'efron'.
#> Start repeat round 3 fold 5 
#> Warning: Starting in glmnet 5.1, the default Cox tie-handling method will change from 'breslow' to 'efron' (matching survival::coxph). To silence this message and lock in the v5.0 default, pass cox.ties = 'breslow' explicitly. To preview the v5.1 behavior, pass cox.ties = 'efron'.

# bootstrap-based discrimination curves has a very narrow band
print(val.boot)
#> High-Dimensional Cox Model Validation Object
#> Random seed: 1010 
#> Validation method: bootstrap
#> Bootstrap samples: 3 
#> Model type: lasso 
#> glmnet model alpha: 1 
#> glmnet model lambda: 0.01839927 
#> glmnet model penalty factor: not specified
#> Time-dependent AUC type: UNO 
#> Evaluation time points for tAUC: 91.25 182.5 273.75 365 456.25 547.5 638.75 730
summary(val.boot)
#> Time-Dependent AUC Summary at Evaluation Time Points
#>              91.25     182.5    273.75       365    456.25     547.5    638.75
#> Mean     0.6039057 0.6833737 0.7393081 0.7802076 0.7804790 0.7930654 0.7652258
#> Min      0.5648485 0.6534571 0.7234113 0.7610369 0.7615473 0.7700398 0.7523180
#> 0.25 Qt. 0.5723232 0.6663138 0.7311501 0.7723847 0.7727783 0.7867470 0.7615195
#> Median   0.5797980 0.6791706 0.7388889 0.7837326 0.7840092 0.8034542 0.7707209
#> 0.75 Qt. 0.6234343 0.6983320 0.7472565 0.7897930 0.7899448 0.8045782 0.7716796
#> Max      0.6670707 0.7174934 0.7556241 0.7958535 0.7958804 0.8057022 0.7726383
#>                730
#> Mean     0.7598040
#> Min      0.7517281
#> 0.25 Qt. 0.7531567
#> Median   0.7545853
#> 0.75 Qt. 0.7638420
#> Max      0.7730988
plot(val.boot)
#>              91.25     182.5    273.75       365    456.25     547.5    638.75
#> Mean     0.6039057 0.6833737 0.7393081 0.7802076 0.7804790 0.7930654 0.7652258
#> Min      0.5648485 0.6534571 0.7234113 0.7610369 0.7615473 0.7700398 0.7523180
#> 0.25 Qt. 0.5723232 0.6663138 0.7311501 0.7723847 0.7727783 0.7867470 0.7615195
#> Median   0.5797980 0.6791706 0.7388889 0.7837326 0.7840092 0.8034542 0.7707209
#> 0.75 Qt. 0.6234343 0.6983320 0.7472565 0.7897930 0.7899448 0.8045782 0.7716796
#> Max      0.6670707 0.7174934 0.7556241 0.7958535 0.7958804 0.8057022 0.7726383
#>                730
#> Mean     0.7598040
#> Min      0.7517281
#> 0.25 Qt. 0.7531567
#> Median   0.7545853
#> 0.75 Qt. 0.7638420
#> Max      0.7730988


# k-fold cv provides a more strict evaluation than bootstrap
print(val.cv)
#> High-Dimensional Cox Model Validation Object
#> Random seed: 1010 
#> Validation method: k-fold cross-validation
#> Cross-validation folds: 5 
#> Model type: lasso 
#> glmnet model alpha: 1 
#> glmnet model lambda: 0.01839927 
#> glmnet model penalty factor: not specified
#> Time-dependent AUC type: UNO 
#> Evaluation time points for tAUC: 91.25 182.5 273.75 365 456.25 547.5 638.75 730
summary(val.cv)
#> Time-Dependent AUC Summary at Evaluation Time Points
#>              91.25     182.5    273.75       365    456.25     547.5    638.75
#> Mean     0.3573288 0.5672998 0.7347406 0.7908407 0.7917342 0.8018140 0.7884141
#> Min      0.0050000 0.0050000 0.5733084 0.5733441 0.5733441 0.5719349 0.5078634
#> 0.25 Qt. 0.0050000 0.4170932 0.6666667 0.7568325 0.7612997 0.7875473 0.7747474
#> Median   0.3571429 0.7057292 0.7263533 0.7837907 0.7837907 0.8128328 0.8227046
#> 0.75 Qt. 0.5306122 0.7705323 0.7705323 0.9033944 0.9033944 0.9012712 0.9012712
#> Max      0.8888889 0.9381443 0.9368421 0.9368421 0.9368421 0.9354839 0.9354839
#>                730
#> Mean     0.7669701
#> Min      0.5078634
#> 0.25 Qt. 0.7731366
#> Median   0.8030970
#> 0.75 Qt. 0.8494822
#> Max      0.9012712
plot(val.cv)
#>              91.25     182.5    273.75       365    456.25     547.5    638.75
#> Mean     0.3573288 0.5672998 0.7347406 0.7908407 0.7917342 0.8018140 0.7884141
#> Min      0.0050000 0.0050000 0.5733084 0.5733441 0.5733441 0.5719349 0.5078634
#> 0.25 Qt. 0.0050000 0.4170932 0.6666667 0.7568325 0.7612997 0.7875473 0.7747474
#> Median   0.3571429 0.7057292 0.7263533 0.7837907 0.7837907 0.8128328 0.8227046
#> 0.75 Qt. 0.5306122 0.7705323 0.7705323 0.9033944 0.9033944 0.9012712 0.9012712
#> Max      0.8888889 0.9381443 0.9368421 0.9368421 0.9368421 0.9354839 0.9354839
#>                730
#> Mean     0.7669701
#> Min      0.5078634
#> 0.25 Qt. 0.7731366
#> Median   0.8030970
#> 0.75 Qt. 0.8494822
#> Max      0.9012712


# repeated cv provides similar results as k-fold cv
# but more robust than k-fold cv
print(val.repcv)
#> High-Dimensional Cox Model Validation Object
#> Random seed: 1010 
#> Validation method: repeated cross-validation
#> Cross-validation folds: 5 
#> Cross-validation repeated times: 3 
#> Model type: lasso 
#> glmnet model alpha: 1 
#> glmnet model lambda: 0.01839927 
#> glmnet model penalty factor: not specified
#> Time-dependent AUC type: UNO 
#> Evaluation time points for tAUC: 91.25 182.5 273.75 365 456.25 547.5 638.75 730
summary(val.repcv)
#> Note: for repeated CV, we evaluated quantile statistic tables for
#> each CV repeat, then calculated element-wise mean across all tables.
#> Time-Dependent AUC Summary at Evaluation Time Points
#>                      91.25     182.5    273.75       365    456.25     547.5
#> Mean of Mean     0.3395427 0.6321354 0.7363690 0.7814654 0.7817636 0.7896420
#> Mean of Min      0.0050000 0.2562963 0.5485496 0.5854828 0.5854828 0.6121502
#> Mean of 0.25 Qt. 0.0050000 0.5178877 0.6777625 0.7360496 0.7372877 0.7383443
#> Mean of Median   0.3707311 0.6772799 0.7465876 0.7788345 0.7788345 0.7803462
#> Mean of 0.75 Qt. 0.5324707 0.8123997 0.8181241 0.8924927 0.8927458 0.8897212
#> Mean of Max      0.7845118 0.8968137 0.8908213 0.9144675 0.9144675 0.9276482
#>                     638.75       730
#> Mean of Mean     0.7652450 0.7615084
#> Mean of Min      0.5378789 0.5378789
#> Mean of 0.25 Qt. 0.6797544 0.6979878
#> Mean of Median   0.7824805 0.7895895
#> Mean of 0.75 Qt. 0.8932617 0.8561390
#> Mean of Max      0.9328494 0.9259470
plot(val.repcv)
#>                      91.25     182.5    273.75       365    456.25     547.5
#> Mean of Mean     0.3395427 0.6321354 0.7363690 0.7814654 0.7817636 0.7896420
#> Mean of Min      0.0050000 0.2562963 0.5485496 0.5854828 0.5854828 0.6121502
#> Mean of 0.25 Qt. 0.0050000 0.5178877 0.6777625 0.7360496 0.7372877 0.7383443
#> Mean of Median   0.3707311 0.6772799 0.7465876 0.7788345 0.7788345 0.7803462
#> Mean of 0.75 Qt. 0.5324707 0.8123997 0.8181241 0.8924927 0.8927458 0.8897212
#> Mean of Max      0.7845118 0.8968137 0.8908213 0.9144675 0.9144675 0.9276482
#>                     638.75       730
#> Mean of Mean     0.7652450 0.7615084
#> Mean of Min      0.5378789 0.5378789
#> Mean of 0.25 Qt. 0.6797544 0.6979878
#> Mean of Median   0.7824805 0.7895895
#> Mean of 0.75 Qt. 0.8932617 0.8561390
#> Mean of Max      0.9328494 0.9259470

# # Test fused lasso, SCAD, and Mnet models
#
# data(smart)
# x = as.matrix(smart[, -c(1, 2)])[1:500,]
# time = smart$TEVENT[1:500]
# event = smart$EVENT[1:500]
# y = survival::Surv(time, event)
#
# set.seed(1010)
# val.boot = validate(
#   x, time, event, model.type = "flasso",
#   lambda1 = 5, lambda2 = 2,
#   method = "bootstrap", boot.times = 10,
#   tauc.type = "UNO", tauc.time = seq(0.25, 2, 0.25) * 365,
#   seed = 1010)
#
# val.cv = validate(
#   x, time, event, model.type = "scad",
#   gamma = 3.7, alpha = 1, lambda = 0.05,
#   method = "cv", nfolds = 5,
#   tauc.type = "UNO", tauc.time = seq(0.25, 2, 0.25) * 365,
#   seed = 1010)
#
# val.repcv = validate(
#   x, time, event, model.type = "mnet",
#   gamma = 3, alpha = 0.3, lambda = 0.05,
#   method = "repeated.cv", nfolds = 5, rep.times = 3,
#   tauc.type = "UNO", tauc.time = seq(0.25, 2, 0.25) * 365,
#   seed = 1010)
#
# print(val.boot)
# summary(val.boot)
# plot(val.boot)
#
# print(val.cv)
# summary(val.cv)
# plot(val.cv)
#
# print(val.repcv)
# summary(val.repcv)
# plot(val.repcv)
```
