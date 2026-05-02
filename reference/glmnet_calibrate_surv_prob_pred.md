# Compute glmnet predicted survival probabilities for calibration

Compute glmnet predicted survival probabilities for calibration

## Usage

``` r
glmnet_calibrate_surv_prob_pred(
  x_tr,
  x_te,
  y_tr,
  alpha,
  lambda,
  pen.factor,
  pred.at
)
```

## Value

list containing predicted survival probability
