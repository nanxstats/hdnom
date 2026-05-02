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
#>  [1]  2  6  8 11 17 21 22 23 24 26
#> 
#> $name
#>  [1] "AGE"      "AAA"      "STENOSIS" "SYSTH"    "HDL"      "GLUT"    
#>  [7] "CREAT"    "IMT"      "ALBUMIN"  "PACKYRS" 
#> 
#> $type
#>  [1] "categorical" "logical"     "logical"     "categorical" "continuous" 
#>  [6] "continuous"  "categorical" "continuous"  "categorical" "continuous" 
#> 
#> $domain
#> $domain[[1]]
#> [1] 19 82
#> 
#> $domain[[2]]
#> [1] 0 1
#> 
#> $domain[[3]]
#> [1] 0 1
#> 
#> $domain[[4]]
#> [1]  79 244
#> 
#> $domain[[5]]
#> [1] 0.22 3.92
#> 
#> $domain[[6]]
#> [1]  2.5 30.9
#> 
#> $domain[[7]]
#> [1]   33 1343
#> 
#> $domain[[8]]
#> [1] 0.36 4.52
#> 
#> $domain[[9]]
#> [1] 1 3
#> 
#> $domain[[10]]
#> [1]   0 120
#> 
#> 
#> attr(,"class")
#> [1] "hdnom.variable.type"
```
