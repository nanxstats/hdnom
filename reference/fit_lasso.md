# Model selection for high-dimensional Cox models with lasso penalty

Automatic model selection for high-dimensional Cox models with lasso
penalty, evaluated by penalized partial-likelihood.

## Usage

``` r
fit_lasso(x, y, nfolds = 5L, rule = c("lambda.min", "lambda.1se"), seed = 1001)
```

## Arguments

- x:

  Data matrix.

- y:

  Response matrix made by
  [`Surv`](https://rdrr.io/pkg/survival/man/Surv.html).

- nfolds:

  Fold numbers of cross-validation.

- rule:

  Model selection criterion, `"lambda.min"` or `"lambda.1se"`. See
  [`cv.glmnet`](https://rdrr.io/pkg/glmnet/man/cv.glmnet.html) for
  details.

- seed:

  A random seed for cross-validation fold division.

## Examples

``` r
data("smart")
x <- as.matrix(smart[, -c(1, 2)])
time <- smart$TEVENT
event <- smart$EVENT
y <- survival::Surv(time, event)

fit <- fit_lasso(x, y, nfolds = 3, rule = "lambda.min", seed = 11)
#> Warning: Starting in glmnet 5.1, the default Cox tie-handling method will change from 'breslow' to 'efron' (matching survival::coxph). To silence this message and lock in the v5.0 default, pass cox.ties = 'breslow' explicitly. To preview the v5.1 behavior, pass cox.ties = 'efron'.
#> Warning: Starting in glmnet 5.1, the default Cox tie-handling method will change from 'breslow' to 'efron' (matching survival::coxph). To silence this message and lock in the v5.0 default, pass cox.ties = 'breslow' explicitly. To preview the v5.1 behavior, pass cox.ties = 'efron'.
#> Warning: Starting in glmnet 5.1, the default Cox tie-handling method will change from 'breslow' to 'efron' (matching survival::coxph). To silence this message and lock in the v5.0 default, pass cox.ties = 'breslow' explicitly. To preview the v5.1 behavior, pass cox.ties = 'efron'.
#> Warning: Starting in glmnet 5.1, the default Cox tie-handling method will change from 'breslow' to 'efron' (matching survival::coxph). To silence this message and lock in the v5.0 default, pass cox.ties = 'breslow' explicitly. To preview the v5.1 behavior, pass cox.ties = 'efron'.
#> Warning: Starting in glmnet 5.1, the default Cox tie-handling method will change from 'breslow' to 'efron' (matching survival::coxph). To silence this message and lock in the v5.0 default, pass cox.ties = 'breslow' explicitly. To preview the v5.1 behavior, pass cox.ties = 'efron'.

nom <- as_nomogram(
  fit, x, time, event,
  pred.at = 365 * 2,
  funlabel = "2-Year Overall Survival Probability"
)

plot(nom)
```
