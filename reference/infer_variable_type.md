# Extract information of selected variables from high-dimensional Cox models

Extract the names and type of selected variables from fitted
high-dimensional Cox models.

## Usage

``` r
infer_variable_type(object, x)
```

## Arguments

- object:

  Model object.

- x:

  Data matrix used to fit the model.

## Value

A list containing the index, name, type and range of the selected
variables.

## Examples

``` r
data("smart")
x <- as.matrix(smart[, -c(1, 2)])
time <- smart$TEVENT
event <- smart$EVENT
y <- survival::Surv(time, event)

fit <- fit_lasso(x, y, nfolds = 3, rule = "lambda.min", seed = 11)
infer_variable_type(fit, x)
#> $index
#>  [1]  1  2  5  6  7  8 11 14 15 17 21 22 23 24 26 27
#> 
#> $name
#>  [1] "SEX"      "AGE"      "CARDIAC"  "AAA"      "PERIPH"   "STENOSIS"
#>  [7] "SYSTH"    "WEIGHT"   "BMI"      "HDL"      "GLUT"     "CREAT"   
#> [13] "IMT"      "ALBUMIN"  "PACKYRS"  "ALCOHOL" 
#> 
#> $type
#>  [1] "logical"     "categorical" "logical"     "logical"     "logical"    
#>  [6] "logical"     "categorical" "categorical" "continuous"  "continuous" 
#> [11] "continuous"  "categorical" "continuous"  "categorical" "continuous" 
#> [16] "categorical"
#> 
#> $domain
#> $domain[[1]]
#> [1] 1 2
#> 
#> $domain[[2]]
#> [1] 19 82
#> 
#> $domain[[3]]
#> [1] 1 0
#> 
#> $domain[[4]]
#> [1] 0 1
#> 
#> $domain[[5]]
#> [1] 1 0
#> 
#> $domain[[6]]
#> [1] 0 1
#> 
#> $domain[[7]]
#> [1]  79 244
#> 
#> $domain[[8]]
#> [1]  37 143
#> 
#> $domain[[9]]
#> [1] 15.43 48.85
#> 
#> $domain[[10]]
#> [1] 0.22 3.92
#> 
#> $domain[[11]]
#> [1]  2.5 30.9
#> 
#> $domain[[12]]
#> [1]   33 1343
#> 
#> $domain[[13]]
#> [1] 0.36 4.52
#> 
#> $domain[[14]]
#> [1] 1 3
#> 
#> $domain[[15]]
#> [1]   0 120
#> 
#> $domain[[16]]
#> [1] 1 3
#> 
#> 
#> attr(,"class")
#> [1] "hdnom.variable.type"
```
