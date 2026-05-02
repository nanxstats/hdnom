# Print high-dimensional Cox model objects

Print high-dimensional Cox model objects

## Usage

``` r
# S3 method for class 'hdnom.model'
print(x, ...)
```

## Arguments

- x:

  Model object.

- ...:

  Other parameters (not used).

## Examples

``` r
data("smart")
x <- as.matrix(smart[, -c(1, 2)])
time <- smart$TEVENT
event <- smart$EVENT
y <- survival::Surv(time, event)

fit <- fit_lasso(x, y, nfolds = 5, rule = "lambda.1se", seed = 11)
print(fit)
#> High-Dimensional Cox Model Object
#> Random seed: 11 
#> Model type: lasso
#> Best lambda: 0.05432858 
```
