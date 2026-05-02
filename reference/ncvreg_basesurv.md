# Breslow baseline hazard estimator for ncvreg objects

Derived from `peperr:::basesurv` and `gbm::basehaz.gbm`.

## Usage

``` r
ncvreg_basesurv(time, event, lp, times.eval = NULL, centered = FALSE)
```

## Arguments

- time:

  Survival time

- event:

  Status indicator

- lp:

  Linear predictors

- times.eval:

  Survival time to evaluate

- centered:

  Should we center the survival curve? See
  [`basehaz`](https://rdrr.io/pkg/survival/man/basehaz.html) for
  details.

## Value

list containing cumulative base hazard

## Examples

``` r
NULL
#> NULL
```
