# Package index

## Building Cox models

Functions for building penalized Cox models.

- [`fit_lasso()`](https://nanx.me/hdnom/reference/fit_lasso.md) : Model
  selection for high-dimensional Cox models with lasso penalty
- [`fit_alasso()`](https://nanx.me/hdnom/reference/fit_alasso.md) :
  Model selection for high-dimensional Cox models with adaptive lasso
  penalty
- [`fit_enet()`](https://nanx.me/hdnom/reference/fit_enet.md) : Model
  selection for high-dimensional Cox models with elastic-net penalty
- [`fit_aenet()`](https://nanx.me/hdnom/reference/fit_aenet.md) : Model
  selection for high-dimensional Cox models with adaptive elastic-net
  penalty
- [`fit_scad()`](https://nanx.me/hdnom/reference/fit_scad.md) : Model
  selection for high-dimensional Cox models with SCAD penalty
- [`fit_snet()`](https://nanx.me/hdnom/reference/fit_snet.md) : Model
  selection for high-dimensional Cox models with Snet penalty
- [`fit_mcp()`](https://nanx.me/hdnom/reference/fit_mcp.md) : Model
  selection for high-dimensional Cox models with MCP penalty
- [`fit_mnet()`](https://nanx.me/hdnom/reference/fit_mnet.md) : Model
  selection for high-dimensional Cox models with Mnet penalty
- [`fit_flasso()`](https://nanx.me/hdnom/reference/fit_flasso.md) :
  Model selection for high-dimensional Cox models with fused lasso
  penalty
- [`print(`*`<hdnom.model>`*`)`](https://nanx.me/hdnom/reference/print.hdnom.model.md)
  : Print high-dimensional Cox model objects
- [`predict(`*`<hdnom.model>`*`)`](https://nanx.me/hdnom/reference/predict.hdnom.model.md)
  : Make predictions from high-dimensional Cox models
- [`infer_variable_type()`](https://nanx.me/hdnom/reference/infer_variable_type.md)
  : Extract information of selected variables from high-dimensional Cox
  models

## Nomogram visualization

Functions for nomogram visualization of the penalized Cox models

- [`as_nomogram()`](https://nanx.me/hdnom/reference/as_nomogram.md) :
  Construct nomogram ojects for high-dimensional Cox models
- [`print(`*`<hdnom.nomogram>`*`)`](https://nanx.me/hdnom/reference/print.hdnom.nomogram.md)
  : Print nomograms objects
- [`plot(`*`<hdnom.nomogram>`*`)`](https://nanx.me/hdnom/reference/plot.hdnom.nomogram.md)
  : Plot nomogram objects

## Model validation

Functions for model validation using bootstrap, k-fold cross-validation,
and repeated k-fold cross-validation. Model performance is assessed by
time-dependent AUC (tAUC).

- [`validate()`](https://nanx.me/hdnom/reference/validate.md) : Validate
  high-dimensional Cox models with time-dependent AUC
- [`print(`*`<hdnom.validate>`*`)`](https://nanx.me/hdnom/reference/print.hdnom.validate.md)
  : Print validation results
- [`summary(`*`<hdnom.validate>`*`)`](https://nanx.me/hdnom/reference/summary.hdnom.validate.md)
  : Summary of validation results
- [`plot(`*`<hdnom.validate>`*`)`](https://nanx.me/hdnom/reference/plot.hdnom.validate.md)
  : Plot optimism-corrected time-dependent discrimination curves for
  validation
- [`validate_external()`](https://nanx.me/hdnom/reference/validate_external.md)
  : Externally validate high-dimensional Cox models with time-dependent
  AUC
- [`print(`*`<hdnom.validate.external>`*`)`](https://nanx.me/hdnom/reference/print.hdnom.validate.external.md)
  : Print external validation results
- [`summary(`*`<hdnom.validate.external>`*`)`](https://nanx.me/hdnom/reference/summary.hdnom.validate.external.md)
  : Summary of external validation results
- [`plot(`*`<hdnom.validate.external>`*`)`](https://nanx.me/hdnom/reference/plot.hdnom.validate.external.md)
  : Plot time-dependent discrimination curves for external validation

## Model calibration

Functions for model calibration using direct fitting, bootstrap
resampling, k-fold cross-validation, and repeated cross-validation.

- [`calibrate()`](https://nanx.me/hdnom/reference/calibrate.md) :
  Calibrate high-dimensional Cox models
- [`print(`*`<hdnom.calibrate>`*`)`](https://nanx.me/hdnom/reference/print.hdnom.calibrate.md)
  : Print calibration results
- [`summary(`*`<hdnom.calibrate>`*`)`](https://nanx.me/hdnom/reference/summary.hdnom.calibrate.md)
  : Summary of calibration results
- [`plot(`*`<hdnom.calibrate>`*`)`](https://nanx.me/hdnom/reference/plot.hdnom.calibrate.md)
  : Plot calibration results
- [`calibrate_external()`](https://nanx.me/hdnom/reference/calibrate_external.md)
  : Externally calibrate high-dimensional Cox models
- [`print(`*`<hdnom.calibrate.external>`*`)`](https://nanx.me/hdnom/reference/print.hdnom.calibrate.external.md)
  : Print external calibration results
- [`summary(`*`<hdnom.calibrate.external>`*`)`](https://nanx.me/hdnom/reference/summary.hdnom.calibrate.external.md)
  : Summary of external calibration results
- [`plot(`*`<hdnom.calibrate.external>`*`)`](https://nanx.me/hdnom/reference/plot.hdnom.calibrate.external.md)
  : Plot external calibration results
- [`kmplot()`](https://nanx.me/hdnom/reference/kmplot.md) : Kaplan-Meier
  plot with number at risk table for internal calibration and external
  calibration results
- [`logrank_test()`](https://nanx.me/hdnom/reference/logrank_test.md) :
  Log-rank test for internal calibration and external calibration
  results

## Model comparison

Functions for model comparison in terms of validation and calibration
performance.

- [`compare_by_validate()`](https://nanx.me/hdnom/reference/compare_by_validate.md)
  : Compare high-dimensional Cox models by model validation
- [`print(`*`<hdnom.compare.validate>`*`)`](https://nanx.me/hdnom/reference/print.hdnom.compare.validate.md)
  : Print model comparison by validation results
- [`summary(`*`<hdnom.compare.validate>`*`)`](https://nanx.me/hdnom/reference/summary.hdnom.compare.validate.md)
  : Summary of model comparison by validation results
- [`plot(`*`<hdnom.compare.validate>`*`)`](https://nanx.me/hdnom/reference/plot.hdnom.compare.validate.md)
  : Plot model comparison by validation results
- [`compare_by_calibrate()`](https://nanx.me/hdnom/reference/compare_by_calibrate.md)
  : Compare high-dimensional Cox models by model calibration
- [`print(`*`<hdnom.compare.calibrate>`*`)`](https://nanx.me/hdnom/reference/print.hdnom.compare.calibrate.md)
  : Print model comparison by calibration results
- [`summary(`*`<hdnom.compare.calibrate>`*`)`](https://nanx.me/hdnom/reference/summary.hdnom.compare.calibrate.md)
  : Summary of model comparison by calibration results
- [`plot(`*`<hdnom.compare.calibrate>`*`)`](https://nanx.me/hdnom/reference/plot.hdnom.compare.calibrate.md)
  : Plot model comparison by calibration results

## Miscellaneous

Miscellaneous functions for supporting survival analysis, such as
baseline hazard estimation and survival curve prediction.

- [`theme_hdnom()`](https://nanx.me/hdnom/reference/theme_hdnom.md) :
  Plot theme (ggplot2) for hdnom
- [`glmnet_basesurv()`](https://nanx.me/hdnom/reference/glmnet_basesurv.md)
  : Breslow baseline hazard estimator for glmnet objects
- [`glmnet_survcurve()`](https://nanx.me/hdnom/reference/glmnet_survcurve.md)
  : Survival curve prediction for glmnet objects
- [`ncvreg_basesurv()`](https://nanx.me/hdnom/reference/ncvreg_basesurv.md)
  : Breslow baseline hazard estimator for ncvreg objects
- [`ncvreg_survcurve()`](https://nanx.me/hdnom/reference/ncvreg_survcurve.md)
  : Survival curve prediction for ncvreg objects
- [`penalized_basesurv()`](https://nanx.me/hdnom/reference/penalized_basesurv.md)
  : Breslow baseline hazard estimator for penfit objects
- [`penalized_survcurve()`](https://nanx.me/hdnom/reference/penalized_survcurve.md)
  : Survival curve prediction for penfit objects
- [`hdnom`](https://nanx.me/hdnom/reference/hdnom-package.md)
  [`hdnom-package`](https://nanx.me/hdnom/reference/hdnom-package.md) :
  hdnom: Benchmarking and Visualization Toolkit for Penalized Cox Models

## Datasets

Example datasets

- [`smart`](https://nanx.me/hdnom/reference/smart.md) : Imputed SMART
  study data
- [`smarto`](https://nanx.me/hdnom/reference/smarto.md) : Original SMART
  study data
