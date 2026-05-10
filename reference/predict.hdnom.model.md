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

fit <- fit_lasso(x, y, nfolds = 5, rule = "lambda.min", seed = 11)
predict(fit, x, y, newx = x[101:105, ], pred.at = 1:10 * 365)
#>            365       730      1095      1460      1825      2190      2555
#> [1,] 0.9168526 0.8745076 0.8260687 0.7800429 0.7228961 0.6772544 0.6220716
#> [2,] 0.9717226 0.9566574 0.9388126 0.9211954 0.8983244 0.8791720 0.8548243
#> [3,] 0.9797314 0.9688642 0.9559284 0.9430883 0.9263137 0.9121723 0.8940669
#> [4,] 0.8633244 0.7969063 0.7236200 0.6566884 0.5773234 0.5169723 0.4476898
#> [5,] 0.9703024 0.9544985 0.9357951 0.9173480 0.8934265 0.8734182 0.8480147
#>           2920      3285      3650
#> [1,] 0.5719206 0.5069896 0.5069896
#> [2,] 0.8314082 0.7989511 0.7989511
#> [3,] 0.8765145 0.8519491 0.8519491
#> [4,] 0.3883074 0.3166438 0.3166438
#> [5,] 0.8236179 0.7898596 0.7898596
```
