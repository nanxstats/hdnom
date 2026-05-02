# Plot model comparison by validation results

Plot model comparison by validation results

## Usage

``` r
# S3 method for class 'hdnom.compare.validate'
plot(
  x,
  interval = FALSE,
  col.pal = c("JCO", "Lancet", "NPG", "AAAS"),
  ylim = NULL,
  ...
)
```

## Arguments

- x:

  An object returned by
  [`compare_by_validate`](https://nanx.me/hdnom/reference/compare_by_validate.md).

- interval:

  Show maximum, minimum, 0.25, and 0.75 quantiles of time-dependent AUC
  as ribbons? Default is `FALSE`.

- col.pal:

  Color palette to use. Possible values are `"JCO"`, `"Lancet"`,
  `"NPG"`, and `"AAAS"`. Default is `"JCO"`.

- ylim:

  Range of y coordinates. For example, `c(0.5, 1)`.

- ...:

  Other parameters (not used).

## Examples

``` r
NULL
#> NULL
```
