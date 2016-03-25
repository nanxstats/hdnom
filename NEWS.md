# CHANGES IN hdnom VERSION 3.6 (2016-03-24)

## NEW FEATURES

* Added 4 new color palettes (JCO, Lancet, NPG, AAAS) for all plots.
  See the vignette for details.

## BUG FIXES

* Fixed no axis problem in Kaplan-Meier plot due to ggplot2 2.1.0 update

## IMPROVEMENTS

* New CSS style for the HTML vignette

# CHANGES IN hdnom VERSION 3.0 (2016-01-03)

## NEW FEATURES

### Package

* New function `hdnom.compare.validate()` for model comparison by validation
* New function `hdnom.compare.calibrate()` for model comparison by calibration
* New function `hdnom.external.validate()` for external validation
* New function `hdnom.external.calibrate()` for external calibration
* New `predict` and `print` methods for `hdcox.model` objects
* New function `hdnom.kmplot()`: Kaplan-Meier analysis for risk groups using
  internal/external calibration results
* New function `hdnom.logrank()`: Log-rank test for risk groups using
  internal/external calibration results

### Web App

* The web app is substantially improved to reflect the new features
* Record random seeds in the generated reports to improve reproducibility
* Allow users to download the model objects for further exploration in `R`

## IMPROVEMENTS

* Improvements on random seed handling for internal validation and calibration
* The vignette is extended to reflect the new features

## BUG FIXES

* Fixed an error in internal calibration which mistakenly used
  testing data when should use training data in computation

# CHANGES IN hdnom VERSION 2.1 (2015-10-26)

## NEW FEATURES

* New website: http://hdnom.org
* Shiny-based web application (http://hdnom.io) shipped.

## IMPROVEMENTS

* Added exception handling for null models in all `hdcox.*()` functions.
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
