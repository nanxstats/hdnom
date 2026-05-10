# Compute validation measures for glmnet objects

Compute validation measures for glmnet objects

## Usage

``` r
glmnet_validate_tauc(
  x_tr,
  x_te,
  y_tr,
  y_te,
  alpha,
  lambda,
  pen.factor,
  tauc.type,
  tauc.time,
  cox.ties
)
```

## Value

time-dependent AUC (tAUC) value
