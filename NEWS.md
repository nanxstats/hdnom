# hdnom 6.0.4

## Bug fixes

- Fixed C code compatibility issues with `_R_USE_STRICT_R_HEADERS_`
  (upcoming default in R 4.5.0) (#13):
  - Defined `STRICT_R_HEADERS`.
  - Replaced `Calloc` and `Free` with `R_Calloc` and `R_Free`.
  - Included `float.h` for `FLT_EPSILON`.

# hdnom 6.0.3

## Improvements

- Fix "lost braces" check notes on r-devel.
- Export S3 method per guidelines in roxygen2 7.3.0.
- Use GitHub Actions to build the pkgdown site.
- Fix code linting issues.

# hdnom 6.0.2

## Improvements

- Use all samples instead of a much smaller set of samples in some code examples
  as a temporary fix for reverse dependency check errors which result in the
  null model. A proper fix will involve updating the logic to use `ncvsurv()`
  in a pathwise manner instead of with a single value of lambda.

# hdnom 6.0.1

## Improvements

- Removed the dependency on `survAUC`. Ported the essential C code for computing
time-dependent AUC and fixed the build issues in r-devel.

# hdnom 6.0.0

This version is a major refactor of the package, with several technical adjustments to improve the functional interface, code structure, and execution performance. As a result, a few critical API-breaking changes have been made. Please update your previous code that calls hdnom accordingly. For the detailed changes, please check the updated items below.

## Improvements

### General

- Renamed exported functions. Most of the exported function have been renamed to be more meaningful and succinct. For example, `hdcox.*()` are renamed as `fit_*()`, `hdnom.nomogram()` is renamed as `as_nomogram()`, `hdnom.validate()` is renamed as `validate()`, and so on.
- Removed the dependency on `rms`, by reusing a minimal set of code from `rms` for nomogram construction and plotting. This results in clearer code structure, better maintainability, and faster package installation/loading speed. Also removed other non-essential package dependencies.
- The first argument for `print` functions are now returned invisbily, to make it easier to use them in a pipe.

### Model Fitting

- The components in the model fitting function returns are now unified across model types. For example, the model object can all be accessed by `fit$model`, and the selected "optimal" hyperparameters can be accessed by `fit$lambda`. The model type is now stored explicitly as `fit$type`.

### Nomograms

- `as_nomogram` (previously `hdnom.nomogram()`) now accepts the fitted model objects directly instead of the `$model` component. It now will recognize the model type automatically, thus the previous arguments `model.type` has been deprecated. so that it is easier to chain the function calls together using `magrittr`.
- In `as_nomogram`, the previous `ddist` argument is not needed anymore and has been removed. There is also no more need to set a `datadist` object as a into the global options variable (which was required in the `rms` user flow).
- The new nomogram implementation prints and plots the nomogram for the penalized regression models directly. This supersedes the old implementation, which fits an OLS model to regress the linear predictors on the same set of predictors selected by the penalized Cox regression model, aiming to approximate the penalized model. The numerical or visual difference is minimal, though.

### Visualizations

- Add a new ggplot2 theme `theme_hdnom()` and applies it to most of the validation, calibration, and comparison plots for a consistent, cleaner look across plots within the package.

# hdnom 5.0 (2018-05-13)

## Improvements

- Exported the survival curve prediction functions (`glmnet.survcurve()`, `ncvreg.survcurve()`, `penalized.survcurve()`) and Breslow baseline hazard estimator functions (`glmnet.basesurv()`, `ncvreg.basesurv()`, `penalized.basesurv()`).
- New URL for the documentation website: https://nanx.me/hdnom/.

# hdnom 4.9 (2017-09-28)

## Improvements

- Use system font stack instead of Google Fonts in vignettes to avoid pandoc SSL issue.

# hdnom 4.8 (2017-03-25)

## Improvements

- Reduced example running time to less than 10s for `hdnom.calibrate()`.

# hdnom 4.7 (2017-03-24)

## Improvements

- Better code indentation style.
- Reduced example running time for fused lasso.
- Updated gallery images in `README.md`.
- HTTPS enabled for the website.

# hdnom 4.6 (2017-01-07)

## Bug Fixes

- Fixed issues in parameter tuning and cross-validation procedures for
fused lasso models ([afc49c9](https://github.com/nanxstats/hdnom/commit/afc49c9ad952edd5881a3e2f14a3503981d213c7)).
The user-visible change is that two parameters `lambda1` and `lambda2`
instead of a single "lambda" are now required to fit, validate, and
calibrate fused lasso models.

## Improvements

- The argument `lambda` in `hdnom.nomogram` is no longer needed and has
been deprecated.
- Allow users to specify `eps` and `max.iter` for MCP and SCAD penalty
related models. Setting the default values to be `1e-4` and `10000`,
which is consistent with ncvreg 3.8-0.

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
- New website design for project website: consistent with the web application.

# hdnom 4.0 (2016-07-05)

## Improvements

- More concrete examples for several functions.
- Introduce argument `ylim` for `plot.hdnom.validate()`,
  `plot.hdnom.external.validate()`, and `plot.hdnom.compare.validate()`
  ([#4](https://github.com/nanxstats/hdnom/issues/4)).

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

- New documentation website.
- New Shiny web application.

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
