# Log-rank test for internal calibration and external calibration results

Log-rank test for internal calibration and external calibration results

## Usage

``` r
logrank_test(object)
```

## Arguments

- object:

  An object returned by
  [`calibrate`](https://nanx.me/hdnom/reference/calibrate.md) or
  [`calibrate_external`](https://nanx.me/hdnom/reference/calibrate_external.md).

## Examples

``` r
data("smart")
# Use the first 1000 samples as training data
# (the data used for internal validation)
x <- as.matrix(smart[, -c(1, 2)])[1:1000, ]
time <- smart$TEVENT[1:1000]
event <- smart$EVENT[1:1000]

# Take the next 1000 samples as external calibration data
# In practice, usually use data collected in other studies
x_new <- as.matrix(smart[, -c(1, 2)])[1001:2000, ]
time_new <- smart$TEVENT[1001:2000]
event_new <- smart$EVENT[1001:2000]

# Fit Cox model with lasso penalty
fit <- fit_lasso(
  x, survival::Surv(time, event),
  nfolds = 5, rule = "lambda.min", seed = 11
)
#> Warning: Starting in glmnet 5.1, the default Cox tie-handling method will change from 'breslow' to 'efron' (matching survival::coxph). To silence this message and lock in the v5.0 default, pass cox.ties = 'breslow' explicitly. To preview the v5.1 behavior, pass cox.ties = 'efron'.
#> Warning: Starting in glmnet 5.1, the default Cox tie-handling method will change from 'breslow' to 'efron' (matching survival::coxph). To silence this message and lock in the v5.0 default, pass cox.ties = 'breslow' explicitly. To preview the v5.1 behavior, pass cox.ties = 'efron'.
#> Warning: Starting in glmnet 5.1, the default Cox tie-handling method will change from 'breslow' to 'efron' (matching survival::coxph). To silence this message and lock in the v5.0 default, pass cox.ties = 'breslow' explicitly. To preview the v5.1 behavior, pass cox.ties = 'efron'.
#> Warning: Starting in glmnet 5.1, the default Cox tie-handling method will change from 'breslow' to 'efron' (matching survival::coxph). To silence this message and lock in the v5.0 default, pass cox.ties = 'breslow' explicitly. To preview the v5.1 behavior, pass cox.ties = 'efron'.
#> Warning: Starting in glmnet 5.1, the default Cox tie-handling method will change from 'breslow' to 'efron' (matching survival::coxph). To silence this message and lock in the v5.0 default, pass cox.ties = 'breslow' explicitly. To preview the v5.1 behavior, pass cox.ties = 'efron'.
#> Warning: Starting in glmnet 5.1, the default Cox tie-handling method will change from 'breslow' to 'efron' (matching survival::coxph). To silence this message and lock in the v5.0 default, pass cox.ties = 'breslow' explicitly. To preview the v5.1 behavior, pass cox.ties = 'efron'.
#> Warning: Starting in glmnet 5.1, the default Cox tie-handling method will change from 'breslow' to 'efron' (matching survival::coxph). To silence this message and lock in the v5.0 default, pass cox.ties = 'breslow' explicitly. To preview the v5.1 behavior, pass cox.ties = 'efron'.

# Internal calibration
cal.int <- calibrate(
  x, time, event,
  model.type = "lasso",
  alpha = 1, lambda = fit$lambda,
  method = "cv", nfolds = 5,
  pred.at = 365 * 9, ngroup = 3
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

logrank_test(cal.int)
#> Call:
#> survdiff(formula = formula("Surv(time, event) ~ grp"))
#> 
#> n=999, 1 observation deleted due to missingness.
#> 
#>         N Observed Expected (O-E)^2/E (O-E)^2/V
#> grp=1 333      127     67.1     53.58      76.5
#> grp=2 333       58     79.1      5.65       8.7
#> grp=3 333       41     79.8     18.86      29.2
#> 
#>  Chisq= 78.4  on 2 degrees of freedom, p= <2e-16 

# External calibration
cal.ext <- calibrate_external(
  fit, x, time, event,
  x_new, time_new, event_new,
  pred.at = 365 * 5, ngroup = 3
)

logrank_test(cal.ext)
#> Call:
#> survdiff(formula = formula("Surv(time, event) ~ grp"))
#> 
#> n=999, 1 observation deleted due to missingness.
#> 
#>         N Observed Expected (O-E)^2/E (O-E)^2/V
#> grp=1 333       86     44.9     37.51      54.2
#> grp=2 333       37     50.2      3.49       5.3
#> grp=3 333       24     51.8     14.94      23.1
#> 
#>  Chisq= 56.1  on 2 degrees of freedom, p= 7e-13 
```
