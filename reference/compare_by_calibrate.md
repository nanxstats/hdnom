# Compare high-dimensional Cox models by model calibration

Compare high-dimensional Cox models by model calibration

## Usage

``` r
compare_by_calibrate(
  x,
  time,
  event,
  model.type = c("lasso", "alasso", "flasso", "enet", "aenet", "mcp", "mnet", "scad",
    "snet"),
  method = c("fitting", "bootstrap", "cv", "repeated.cv"),
  boot.times = NULL,
  nfolds = NULL,
  rep.times = NULL,
  pred.at,
  ngroup = 5,
  rule = c("lambda.min", "lambda.1se"),
  cox.ties = c("breslow", "efron"),
  seed = 1001,
  trace = TRUE
)
```

## Arguments

- x:

  Matrix of training data used for fitting the model; on which to run
  the calibration.

- time:

  Survival time. Must be of the same length with the number of rows as
  `x`.

- event:

  Status indicator, normally 0 = alive, 1 = dead. Must be of the same
  length with the number of rows as `x`.

- model.type:

  Model types to compare. Could be at least two of `"lasso"`,
  `"alasso"`, `"flasso"`, `"enet"`, `"aenet"`, `"mcp"`, `"mnet"`,
  `"scad"`, or `"snet"`.

- method:

  Calibration method. Could be `"bootstrap"`, `"cv"`, or
  `"repeated.cv"`.

- boot.times:

  Number of repetitions for bootstrap.

- nfolds:

  Number of folds for cross-validation and repeated cross-validation.

- rep.times:

  Number of repeated times for repeated cross-validation.

- pred.at:

  Time point at which calibration should take place.

- ngroup:

  Number of groups to be formed for calibration.

- rule:

  Model selection criterion for glmnet models, \`"lambda.min"\` or
  \`"lambda.1se"\`. Defaults to \`"lambda.min"\`.

- cox.ties:

  Cox tie-handling method for glmnet model fits and refits.

- seed:

  A random seed for cross-validation fold division.

- trace:

  Logical. Output the calibration progress or not. Default is `TRUE`.

## Examples

``` r
data(smart)
x <- as.matrix(smart[, -c(1, 2)])
time <- smart$TEVENT
event <- smart$EVENT

# Compare lasso and adaptive lasso by 5-fold cross-validation
cmp.cal.cv <- compare_by_calibrate(
  x, time, event,
  model.type = c("lasso", "alasso"),
  method = "fitting",
  pred.at = 365 * 9, ngroup = 5, seed = 1001
)
#> Starting model 1 : lasso 
#> Start fitting ...
#> Starting model 2 : alasso 
#> Start fitting ...

print(cmp.cal.cv)
#> High-Dimensional Cox Model Calibration Object
#> Random seed: 1001 
#> Calibration method: fitting
#> Model type: lasso 
#> glmnet model alpha: 1 
#> glmnet model lambda: 0.002521583 
#> glmnet model penalty factor: not specified
#> Calibration time point: 3285 
#> Number of groups formed for calibration: 5 
#> 
#> High-Dimensional Cox Model Calibration Object
#> Random seed: 1001 
#> Calibration method: fitting
#> Model type: alasso 
#> glmnet model alpha: 1 
#> glmnet model lambda: 0.001912884 
#> glmnet model penalty factor: specified
#> Calibration time point: 3285 
#> Number of groups formed for calibration: 5 
#> 
summary(cmp.cal.cv)
#>   Model type: lasso 
#>   Calibration Summary Table
#>   Predicted  Observed Lower 95% Upper 95%
#> 1 0.5466818 0.4637830 0.3857102 0.5576587
#> 2 0.7196071 0.7316329 0.6710281 0.7977112
#> 3 0.7914195 0.8199814 0.7511298 0.8951443
#> 4 0.8386738 0.8604124 0.7892224 0.9380240
#> 5 0.8858013 0.9279530 0.8968975 0.9600839
#> attr(,"cox.ties")
#> [1] "breslow"
#> 
#>   Model type: alasso 
#>   Calibration Summary Table
#>   Predicted  Observed Lower 95% Upper 95%
#> 1 0.5341084 0.4535579 0.3724883 0.5522718
#> 2 0.7251070 0.7419299 0.6812376 0.8080294
#> 3 0.8014661 0.8077740 0.7324507 0.8908434
#> 4 0.8520771 0.8829392 0.8388129 0.9293868
#> 5 0.8996329 0.9265275 0.8950859 0.9590736
#> attr(,"cox.ties")
#> [1] "breslow"
#> 
plot(cmp.cal.cv)
```
