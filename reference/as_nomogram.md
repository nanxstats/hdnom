# Construct nomogram ojects for high-dimensional Cox models

Construct nomograms ojects for high-dimensional Cox models

## Usage

``` r
as_nomogram(
  object,
  x,
  time,
  event,
  pred.at = NULL,
  fun.at = NULL,
  funlabel = NULL
)
```

## Arguments

- object:

  Model object fitted by \`hdnom::fit\_\*()\` functions.

- x:

  Matrix of training data used for fitting the model.

- time:

  Survival time. Must be of the same length with the number of rows as
  `x`.

- event:

  Status indicator, normally 0 = alive, 1 = dead. Must be of the same
  length with the number of rows as `x`.

- pred.at:

  Time point at which to plot nomogram prediction axis.

- fun.at:

  Function values to label on axis.

- funlabel:

  Label for `fun` axis.

## Note

The nomogram visualizes the model under the automatically selected
"optimal" hyperparameters (e.g. lambda, alpha, gamma).

## Examples

``` r
data(smart)
x <- as.matrix(smart[, -c(1, 2)])
time <- smart$TEVENT
event <- smart$EVENT
y <- survival::Surv(time, event)

fit <- fit_lasso(x, y, nfolds = 5, rule = "lambda.1se", seed = 1001)

nom <- as_nomogram(
  fit, x, time, event,
  pred.at = 365 * 2,
  funlabel = "2-Year Overall Survival Probability"
)

print(nom)
#> Points per unit of linear predictor: 107.0109 
#> Linear predictor units per point   : 0.009344843 
#> 
#> 
#>  AGE Points
#>  15   0    
#>  20   4    
#>  25   8    
#>  30  12    
#>  35  16    
#>  40  21    
#>  45  25    
#>  50  29    
#>  55  33    
#>  60  37    
#>  65  41    
#>  70  45    
#>  75  49    
#>  80  53    
#>  85  58    
#> 
#> 
#>  AAA Points
#>  0    0    
#>  1   16    
#> 
#> 
#>  CREAT Points
#>     0    0   
#>   100    7   
#>   200   14   
#>   300   21   
#>   400   29   
#>   500   36   
#>   600   43   
#>   700   50   
#>   800   57   
#>   900   64   
#>  1000   71   
#>  1100   79   
#>  1200   86   
#>  1300   93   
#>  1400  100   
#> 
#> 
#>  IMT Points
#>  0.0  0    
#>  0.5  4    
#>  1.0  7    
#>  1.5 11    
#>  2.0 15    
#>  2.5 18    
#>  3.0 22    
#>  3.5 26    
#>  4.0 29    
#>  4.5 33    
#>  5.0 37    
#> 
#> 
#>  ALBUMIN Points
#>  1        0    
#>  2       18    
#>  3       36    
#> 
#> 
#>  Total Points 2-Year Overall Survival Probability
#>           113                                0.90
#>            36                                0.95
#> 
plot(nom)
```
