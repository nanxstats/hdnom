# hdnom 4.6 (2017-01-05)

## Bug Fixes

- Fixed issues in parameter tuning and cross-validation procedures for
fused lasso models ([afc49c9](https://github.com/road2stat/hdnom/commit/afc49c9ad952edd5881a3e2f14a3503981d213c7)).
The user visible change is that two parameters `lambda1` and `lambda2`
instead of a single "lambda" are now required to fit, validate, and
calibrate fused lasso models.

## Improvements

- The argument `lambda` in `hdnom.nomogram` is no longer needed and has
been deprecated.

# hdnom 4.5 (2016-12-24)

## Bug Fixes

- Fixed vanishing axis problem in Kaplan-Meier plot `hdnom.kmplot()` under
ggplot2 2.2.0, which is caused by a previous workaround for a bug introduced
in ggplot2 2.1.0.
- Fixed potential convergence issues for examples under ncvreg >= 3.7-0 new
convergence criterion, by increasing `max.iter` for `ncvsurv` to a
substantially higher value (5e+4).
- Fixed single lambda support issues in `ncvsurv` under ncvreg >= 3.7-0.

## Improvements

- Added Windows continuous integration using AppVeyor.
- New website design for hdnom.org: consistent with the web application hdnom.io.

# hdnom 4.0 (2016-07-05)

## Improvements

- More concrete examples for several functions.
- Introduce argument `ylim` for `plot.hdnom.validate()`,
  `plot.hdnom.external.validate()`, and `plot.hdnom.compare.validate()`
  ([#4](https://github.com/road2stat/hdnom/issues/4)).

# hdnom 3.7 (2016-03-25)

## Bug Fixes

- Removed one redundant color from the lancet color palette.

# hdnom 3.6 (2016-03-24)

## New Features

- Added 4 new color palettes (JCO, Lancet, NPG, AAAS) for all plots.
  See the vignette for details.

## Bug Fixes

- Fixed vanishing axis problem in Kaplan-Meier plot due to ggplot2 2.1.0 update.

## Improvements

- New CSS style for the HTML vignette.

# hdnom 3.0 (2016-01-03)

## New Features

### Package

- New function `hdnom.compare.validate()` for model comparison by validation
- New function `hdnom.compare.calibrate()` for model comparison by calibration
- New function `hdnom.external.validate()` for external validation
- New function `hdnom.external.calibrate()` for external calibration
- New `predict` and `print` methods for `hdcox.model` objects
- New function `hdnom.kmplot()`: Kaplan-Meier analysis for risk groups using
  internal/external calibration results
- New function `hdnom.logrank()`: Log-rank test for risk groups using
  internal/external calibration results

### Web Application

- The web application is substantially improved to reflect the new features
- Record random seeds in the generated reports to improve reproducibility
- Allow users to download the model objects for further exploration in `R`

## Improvements

- Improvements on random seed handling for internal validation and calibration
- The vignette is extended to reflect the new features

## Bug Fixes

- Fixed an error in internal calibration which mistakenly used
  testing data when should use training data in computation

# hdnom 2.1 (2015-10-26)

## New Features

- New website: http://hdnom.org
- Shiny-based web application (http://hdnom.io) shipped.

## Improvements

- Added exception handling for null models in all `hdcox.*()` functions.
  Make examples compatible with ncvreg 3.5-0, which refined CV
  implementation for survival models and improved computation speed.

# hdnom 2.0 (2015-09-15)

## New Features

- Support five more high-dimensional penalized Cox model types:

  - Fused lasso
  - MCP
  - Mnet
  - SCAD
  - Snet

# hdnom 1.2 (2015-08-27)

## Improvements

- Reduced example running time for `hdnom.validate()`, `hdnom.calibrate()`,
`hdcox.aenet()`, and `hdcox.enet()` by reducing resampling times.

# hdnom 1.1 (2015-08-27)

## Bug Fixes

- Added argument `parallel` to `hdcox.aenet()` and `hdcox.enet()` to enable
or disable the use of parallel parameter tuning.

# hdnom 1.0 (2015-08-04)

## New Features

- Initial version
- Nomograms for glmnet models
- Validation for glmnet models
- Calibration for glmnet models
