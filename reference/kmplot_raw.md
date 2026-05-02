# Kaplan-Meier Plot with Number at Risk Table

Kaplan-Meier Plot with Number at Risk Table

## Usage

``` r
kmplot_raw(
  fit,
  group.name = NULL,
  time.at = NULL,
  surv.df = NULL,
  col.pal = NULL
)
```

## Arguments

- fit:

  a [`survfit`](https://rdrr.io/pkg/survival/man/survfit.html) object

- group.name:

  Group labels. Default is Group 1, Group 2, ... Group k.

- time.at:

  Time points to evaluate the number at risk.

- surv.df:

  Data frame containing survival time, event and risk group for log-rank
  test.

- col.pal:

  color palette to use
