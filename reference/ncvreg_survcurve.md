# Survival curve prediction for ncvreg objects

Derived from c060::predictProb.coxnet

## Usage

``` r
ncvreg_survcurve(object, time, event, x, survtime)
```

## Arguments

- object:

  `ncvreg` model object

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
