# Calibrate high-dimensional Cox models

Calibrate high-dimensional Cox models

## Usage

``` r
calibrate(
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
  method = c("fitting", "bootstrap", "cv", "repeated.cv"),
  boot.times = NULL,
  nfolds = NULL,
  rep.times = NULL,
  pred.at,
  ngroup = 5,
  seed = 1001,
  trace = TRUE
)
```

## Arguments

- x:

  Matrix of training data used for fitting the model; on which to run
  the calibration.

- time:

  Survival time. Must be of the same length with the number of rows as
  `x`.

- event:

  Status indicator, normally 0 = alive, 1 = dead. Must be of the same
  length with the number of rows as `x`.

- model.type:

  Model type to calibrate. Could be one of `"lasso"`, `"alasso"`,
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
  resampled data. From the Cox model you have built.

- pen.factor:

  Penalty factors to apply to each coefficient. From the built *adaptive
  lasso* or *adaptive elastic-net* model.

- gamma:

  Value of the model parameter gamma for MCP/SCAD/Mnet/Snet models.

- lambda1:

  Value of the penalty parameter lambda1 for fused lasso model.

- lambda2:

  Value of the penalty parameter lambda2 for fused lasso model.

- method:

  Calibration method. Options including `"fitting"`, `"bootstrap"`,
  `"cv"`, and `"repeated.cv"`.

- boot.times:

  Number of repetitions for bootstrap.

- nfolds:

  Number of folds for cross-validation and repeated cross-validation.

- rep.times:

  Number of repeated times for repeated cross-validation.

- pred.at:

  Time point at which calibration should take place.

- ngroup:

  Number of groups to be formed for calibration.

- seed:

  A random seed for resampling.

- trace:

  Logical. Output the calibration progress or not. Default is `TRUE`.

## Examples

``` r
data("smart")
x <- as.matrix(smart[, -c(1, 2)])
time <- smart$TEVENT
event <- smart$EVENT
y <- survival::Surv(time, event)

# Fit Cox model with lasso penalty
fit <- fit_lasso(x, y, nfolds = 5, rule = "lambda.1se", seed = 1001)

# Model calibration by fitting the original data directly
cal.fitting <- calibrate(
  x, time, event,
  model.type = "lasso",
  alpha = 1, lambda = fit$lambda,
  method = "fitting",
  pred.at = 365 * 9, ngroup = 5,
  seed = 1010
)
#> Start fitting ...

# Model calibration by 5-fold cross-validation
cal.cv <- calibrate(
  x, time, event,
  model.type = "lasso",
  alpha = 1, lambda = fit$lambda,
  method = "cv", nfolds = 5,
  pred.at = 365 * 9, ngroup = 5,
  seed = 1010
)
#> Start fold 1 
#> Start fold 2 
#> Start fold 3 
#> Start fold 4 
#> Start fold 5 

print(cal.fitting)
#> High-Dimensional Cox Model Calibration Object
#> Random seed: 1010 
#> Calibration method: fitting
#> Model type: lasso 
#> glmnet model alpha: 1 
#> glmnet model lambda: 0.03744658 
#> glmnet model penalty factor: not specified
#> Calibration time point: 3285 
#> Number of groups formed for calibration: 5 
summary(cal.fitting)
#>   Calibration Summary Table
#>   Predicted  Observed Lower 95% Upper 95%
#> 1 0.6609428 0.4601632 0.3836063 0.5519988
#> 2 0.7078401 0.7520519 0.6977012 0.8106366
#> 3 0.7317541 0.7441680 0.6364450 0.8701239
#> 4 0.7512392 0.8413863 0.7660086 0.9241814
#> 5 0.7741434 0.8938572 0.8568490 0.9324639
plot(cal.fitting)
#> Warning: `aes_string()` was deprecated in ggplot2 3.0.0.
#> ℹ Please use tidy evaluation idioms with `aes()`.
#> ℹ See also `vignette("ggplot2-in-packages")` for more information.
#> ℹ The deprecated feature was likely used in the hdnom package.
#>   Please report the issue at <https://github.com/nanxstats/hdnom/issues>.


print(cal.cv)
#> High-Dimensional Cox Model Calibration Object
#> Random seed: 1010 
#> Calibration method: k-fold cross-validation
#> Cross-validation folds: 5 
#> Model type: lasso 
#> glmnet model alpha: 1 
#> glmnet model lambda: 0.03744658 
#> glmnet model penalty factor: not specified
#> Calibration time point: 3285 
#> Number of groups formed for calibration: 5 
summary(cal.cv)
#>   Calibration Summary Table
#>   Predicted  Observed Lower 95% Upper 95%
#> 1 0.6588508 0.5078961 0.4309993 0.5985125
#> 2 0.7099409 0.7212081 0.6604775 0.7875229
#> 3 0.7341094 0.7168790 0.5481497 0.9375459
#> 4 0.7530252 0.8275354 0.7520007 0.9106571
#> 5 0.7745243 0.8569103 0.7987051 0.9193573
plot(cal.cv)


# # Test fused lasso, SCAD, and Mnet models
# data(smart)
# x = as.matrix(smart[, -c(1, 2)])[1:500, ]
# time = smart$TEVENT[1:500]
# event = smart$EVENT[1:500]
# y = survival::Surv(time, event)
#
# set.seed(1010)
# cal.fitting = calibrate(
#   x, time, event, model.type = "flasso",
#   lambda1 = 5, lambda2 = 2,
#   method = "fitting",
#   pred.at = 365 * 9, ngroup = 5,
#   seed = 1010)
#
# cal.boot = calibrate(
#   x, time, event, model.type = "scad",
#   gamma = 3.7, alpha = 1, lambda = 0.03,
#   method = "bootstrap", boot.times = 10,
#   pred.at = 365 * 9, ngroup = 5,
#   seed = 1010)
#
# cal.cv = calibrate(
#   x, time, event, model.type = "mnet",
#   gamma = 3, alpha = 0.3, lambda = 0.03,
#   method = "cv", nfolds = 5,
#   pred.at = 365 * 9, ngroup = 5,
#   seed = 1010)
#
# cal.repcv = calibrate(
#   x, time, event, model.type = "flasso",
#   lambda1 = 5, lambda2 = 2,
#   method = "repeated.cv", nfolds = 5, rep.times = 3,
#   pred.at = 365 * 9, ngroup = 5,
#   seed = 1010)
#
# print(cal.fitting)
# summary(cal.fitting)
# plot(cal.fitting)
#
# print(cal.boot)
# summary(cal.boot)
# plot(cal.boot)
#
# print(cal.cv)
# summary(cal.cv)
# plot(cal.cv)
#
# print(cal.repcv)
# summary(cal.repcv)
# plot(cal.repcv)
```
