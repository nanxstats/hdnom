# CHANGES IN hdnom VERSION 2.1 (2015-10-26)

## NEW FEATURES

* New website: hdnom.org
* Shiny-based web application (hdnom.io) shipped.

## IMPROVEMENTS

* Added exception handling for null models in all hdcox.* functions.
  Make examples compatible with ncvreg 3.5-0, which refined CV
  implementation for survival models and improved computation speed.

# CHANGES IN hdnom VERSION 2.0 (2015-09-15)

## NEW FEATURES

* Support five more high-dimensional penalized Cox model types:

  * Fused lasso
  * MCP
  * Mnet
  * SCAD
  * Snet

# CHANGES IN hdnom VERSION 1.2 (2015-08-27)

## IMPROVEMENTS

* Reduced example running time for `hdnom.validate()`, `hdnom.calibrate()`,
`hdcox.aenet()`, and `hdcox.enet()` by reducing resampling times.

# CHANGES IN hdnom VERSION 1.1 (2015-08-27)

## BUG FIXES

* Added argument `parallel` to `hdcox.aenet()` and `hdcox.enet()` to enable
or disable the use of parallel parameter tuning.

# CHANGES IN hdnom VERSION 1.0 (2015-08-04)

## NEW FEATURES

* Initial version of hdnom.
* Nomograms for glmnet models.
* Validation for glmnet models.
* Calibration for glmnet models.
