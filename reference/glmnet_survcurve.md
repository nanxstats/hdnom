# Survival curve prediction for glmnet objects

Derived from c060::predictProb.coxnet

## Usage

``` r
glmnet_survcurve(object, time, event, x, survtime)
```

## Arguments

- object:

  `glmnet` model object

- time:

  Survival time

- event:

  Status indicator

- x:

  Predictor matrix

- survtime:

  Survival time to evaluate

## Value

list containing predicted survival probabilities and linear predictors
for all samples

## Examples

``` r
NULL
#> NULL
```
