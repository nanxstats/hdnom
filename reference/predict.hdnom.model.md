# Make predictions from high-dimensional Cox models

Predict overall survival probability at certain time points from fitted
Cox models.

## Usage

``` r
# S3 method for class 'hdnom.model'
predict(object, x, y, newx, pred.at, ...)
```

## Arguments

- object:

  Model object.

- x:

  Data matrix used to fit the model.

- y:

  Response matrix made with
  [`Surv`](https://rdrr.io/pkg/survival/man/Surv.html).

- newx:

  Matrix (with named columns) of new values for `x` at which predictions
  are to be made.

- pred.at:

  Time point at which prediction should take place.

- ...:

  Other parameters (not used).

## Value

A `nrow(newx) x length(pred.at)` matrix containing overall survival
probablity.

## Examples

``` r
data("smart")
x <- as.matrix(smart[, -c(1, 2)])
time <- smart$TEVENT
event <- smart$EVENT
y <- survival::Surv(time, event)

fit <- fit_lasso(x, y, nfolds = 5, rule = "lambda.1se", seed = 11)
predict(fit, x, y, newx = x[101:105, ], pred.at = 1:10 * 365)
#>            365       730      1095      1460      1825      2190      2555
#> [1,] 0.9597234 0.9384128 0.9134911 0.8894236 0.8586872 0.8338478 0.8020999
#> [2,] 0.9597234 0.9384128 0.9134911 0.8894236 0.8586872 0.8338478 0.8020999
#> [3,] 0.9597234 0.9384128 0.9134911 0.8894236 0.8586872 0.8338478 0.8020999
#> [4,] 0.9597234 0.9384128 0.9134911 0.8894236 0.8586872 0.8338478 0.8020999
#> [5,] 0.9597234 0.9384128 0.9134911 0.8894236 0.8586872 0.8338478 0.8020999
#>           2920      3285      3650
#> [1,] 0.7705992 0.7254538 0.7254538
#> [2,] 0.7705992 0.7254538 0.7254538
#> [3,] 0.7705992 0.7254538 0.7254538
#> [4,] 0.7705992 0.7254538 0.7254538
#> [5,] 0.7705992 0.7254538 0.7254538
```
