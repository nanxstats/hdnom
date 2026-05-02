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

# Internal calibration
cal.int <- calibrate(
  x, time, event,
  model.type = "lasso",
  alpha = 1, lambda = fit$lambda,
  method = "cv", nfolds = 5,
  pred.at = 365 * 9, ngroup = 3
)
#> Start fold 1 
#> Start fold 2 
#> Start fold 3 
#> Start fold 4 
#> Start fold 5 

logrank_test(cal.int)
#> Call:
#> survdiff(formula = formula("Surv(time, event) ~ grp"))
#> 
#> n=999, 1 observation deleted due to missingness.
#> 
#>         N Observed Expected (O-E)^2/E (O-E)^2/V
#> grp=1 333      127     66.8     54.26     77.31
#> grp=2 333       59     79.2      5.14      7.92
#> grp=3 333       40     80.0     20.02     31.06
#> 
#>  Chisq= 79.7  on 2 degrees of freedom, p= <2e-16 

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
#> grp=1 333       90     44.3     47.14      67.7
#> grp=2 333       32     51.1      7.11      10.9
#> grp=3 333       25     51.6     13.75      21.2
#> 
#>  Chisq= 68.2  on 2 degrees of freedom, p= 2e-15 
```
