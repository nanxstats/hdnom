#' Compare High-Dimensional Cox Models by Model Calibration
#'
#' Compare High-Dimensional Cox Models by Model Calibration
#'
#' @param x Matrix of training data used for fitting the model;
#' on which to run the calibration.
#' @param time Survival time.
#' Must be of the same length with the number of rows as \code{x}.
#' @param event Status indicator, normally 0 = alive, 1 = dead.
#' Must be of the same length with the number of rows as \code{x}.
#' @param model.type Model types to compare. Could be at least two of
#' \code{"lasso"}, \code{"alasso"}, \code{"flasso"}, \code{"enet"},
#' \code{"aenet"}, \code{"mcp"}, \code{"mnet"}, \code{"scad"},
#' or \code{"snet"}.
#' @param method Calibration method.
#' Could be \code{"bootstrap"}, \code{"cv"}, or \code{"repeated.cv"}.
#' @param boot.times Number of repetitions for bootstrap.
#' @param nfolds Number of folds for cross-validation and
#' repeated cross-validation.
#' @param rep.times Number of repeated times for repeated cross-validation.
#' @param pred.at Time point at which calibration should take place.
#' @param ngroup Number of groups to be formed for calibration.
#' @param seed A random seed for cross-validation fold division.
#' @param trace Logical. Output the calibration progress or not.
#' Default is \code{TRUE}.
#'
#' @export hdnom.compare.calibrate
#'
#' @examples
#' # Load imputed SMART data
#' data(smart)
#' x = as.matrix(smart[, -c(1, 2)])
#' time = smart$TEVENT
#' event = smart$EVENT
#'
#' # Compare lasso and adaptive lasso by 5-fold cross-validation
#' cmp.cal.cv =
#'   hdnom.compare.calibrate(x, time, event,
#'                           model.type = c("lasso", "alasso"),
#'                           method = "fitting",
#'                           pred.at = 365 * 9, ngroup = 5, seed = 1001)
#'
#' print(cmp.cal.cv)
#' summary(cmp.cal.cv)
#' plot(cmp.cal.cv)
hdnom.compare.calibrate =
  function(x, time, event,
           model.type = c('lasso', 'alasso', 'flasso',
                          'enet', 'aenet',
                          'mcp', 'mnet',
                          'scad', 'snet'),
           method = c('fitting', 'bootstrap', 'cv', 'repeated.cv'),
           boot.times = NULL, nfolds = NULL, rep.times = NULL,
           pred.at, ngroup = 5,
           seed = 1001, trace = TRUE) {

    method = match.arg(method)
    if (length(pred.at) != 1L) stop('pred.at should only contain 1 time point')

    if (!all(model.type %in% c('lasso', 'alasso', 'flasso', 'enet', 'aenet',
                               'mcp', 'mnet', 'scad', 'snet'))) {
      stop('Unknown model type(s) specified')
    }

    nmodel = length(model.type)
    problist = vector('list', nmodel)

    # check parameters for different methods
    if (method == 'fitting') {
      if (!is.null(boot.times) || !is.null(nfolds) || !is.null(rep.times))
        stop('boot.times, nfolds, and rep.times must be NULL when method = "fitting"')
    }

    if (method == 'bootstrap') {
      if (!is.null(nfolds) || !is.null(rep.times))
        stop('nfolds and rep.times must be NULL when method = "bootstrap"')
      if (is.null(boot.times)) stop('please specify boot.times')
    }

    if (method == 'cv') {
      if (!is.null(boot.times) || !is.null(rep.times))
        stop('boot.times and rep.times must be NULL when method = "cv"')
      if (is.null(nfolds)) stop('please specify nfolds')
    }

    if (method == 'repeated.cv') {
      if (!is.null(boot.times))
        stop('boot.times must be NULL when method = "repeated.cv"')
      if (is.null(nfolds) || is.null(rep.times))
        stop('please specify nfolds and rep.times')
    }

    for (i in 1L:nmodel) {

      if (trace) cat('Starting model', i, ':', model.type[i], '\n')

      switch(model.type[i],

             lasso = {

               cvfit = hdcox.lasso(x, Surv(time, event), nfolds = 5L,
                                   rule = 'lambda.1se', seed = seed)

               problist[[i]] =
                 hdnom.calibrate(
                   x, time, event, model.type = 'lasso',
                   alpha = 1, lambda = cvfit$'lasso_best_lambda',
                   method = method, boot.times = boot.times, nfolds = nfolds, rep.times = rep.times,
                   pred.at = pred.at, ngroup = ngroup,
                   seed = seed, trace = trace)

             },

             alasso = {

               cvfit = hdcox.alasso(x, Surv(time, event), nfolds = 5L,
                                    rule = 'lambda.1se', seed = rep(seed, 2))

               problist[[i]] =
                 hdnom.calibrate(
                   x, time, event, model.type = 'alasso',
                   alpha = 1, lambda = cvfit$'alasso_best_lambda', pen.factor = cvfit$'pen_factor',
                   method = method, boot.times = boot.times, nfolds = nfolds, rep.times = rep.times,
                   pred.at = pred.at, ngroup = ngroup,
                   seed = seed, trace = trace)

             },

             flasso = {

               cvfit = hdcox.flasso(x, Surv(time, event), nfolds = 5L, seed = seed)

               problist[[i]] =
                 hdnom.calibrate(
                   x, time, event, model.type = 'flasso',
                   alpha = 1, lambda = cvfit$'flasso_best_lambda',
                   method = method, boot.times = boot.times, nfolds = nfolds, rep.times = rep.times,
                   pred.at = pred.at, ngroup = ngroup,
                   seed = seed, trace = trace)

             },

             enet = {

               cvfit = hdcox.enet(
                 x, Surv(time, event), nfolds = 5L,
                 alphas = c(0.1, 0.25, 0.5, 0.75, 0.9),  # to reduce computation time
                 rule = 'lambda.1se', seed = seed)

               problist[[i]] =
                 hdnom.calibrate(
                   x, time, event, model.type = 'enet',
                   alpha = cvfit$'enet_best_alpha', lambda = cvfit$'enet_best_lambda',
                   method = method, boot.times = boot.times, nfolds = nfolds, rep.times = rep.times,
                   pred.at = pred.at, ngroup = ngroup,
                   seed = seed, trace = trace)

             },

             aenet = {

               cvfit = hdcox.aenet(
                 x, Surv(time, event), nfolds = 5L,
                 alphas = c(0.1, 0.25, 0.5, 0.75, 0.9),  # to reduce computation time
                 rule = 'lambda.1se', seed = rep(seed, 2))

               problist[[i]] =
                 hdnom.calibrate(
                   x, time, event, model.type = 'aenet',
                   alpha = cvfit$'aenet_best_alpha', lambda = cvfit$'aenet_best_lambda', pen.factor = cvfit$'pen_factor',
                   method = method, boot.times = boot.times, nfolds = nfolds, rep.times = rep.times,
                   pred.at = pred.at, ngroup = ngroup,
                   seed = seed, trace = trace)

             },

             mcp = {

               cvfit = hdcox.mcp(x, Surv(time, event), nfolds = 5L, seed = seed)

               problist[[i]] =
                 hdnom.calibrate(
                   x, time, event, model.type = 'mcp',
                   alpha = 1, gamma = cvfit$'mcp_best_gamma', lambda = cvfit$'mcp_best_lambda',
                   method = method, boot.times = boot.times, nfolds = nfolds, rep.times = rep.times,
                   pred.at = pred.at, ngroup = ngroup,
                   seed = seed, trace = trace)

             },

             mnet = {

               cvfit = hdcox.mnet(
                 x, Surv(time, event), nfolds = 5L,
                 alphas = c(0.1, 0.25, 0.5, 0.75, 0.9), # to reduce computation time
                 seed = seed)

               problist[[i]] =
                 hdnom.calibrate(
                   x, time, event, model.type = 'mnet',
                   alpha = cvfit$'mnet_best_alpha', gamma = cvfit$'mnet_best_gamma', lambda = cvfit$'mnet_best_lambda',
                   method = method, boot.times = boot.times, nfolds = nfolds, rep.times = rep.times,
                   pred.at = pred.at, ngroup = ngroup,
                   seed = seed, trace = trace)

             },

             scad = {

               cvfit = hdcox.scad(x, Surv(time, event), nfolds = 5L, seed = seed)

               problist[[i]] =
                 hdnom.calibrate(
                   x, time, event, model.type = 'scad',
                   alpha = 1, gamma = cvfit$'scad_best_gamma', lambda = cvfit$'scad_best_lambda',
                   method = method, boot.times = boot.times, nfolds = nfolds, rep.times = rep.times,
                   pred.at = pred.at, ngroup = ngroup,
                   seed = seed, trace = trace)

             },

             snet = {

               cvfit = hdcox.snet(
                 x, Surv(time, event), nfolds = 5L,
                 alphas = c(0.1, 0.25, 0.5, 0.75, 0.9), # to reduce computation time
                 seed = seed)

               problist[[i]] =
                 hdnom.calibrate(
                   x, time, event, model.type = 'snet',
                   alpha = cvfit$'snet_best_alpha', gamma = cvfit$'snet_best_gamma', lambda = cvfit$'snet_best_lambda',
                   method = method, boot.times = boot.times, nfolds = nfolds, rep.times = rep.times,
                   pred.at = pred.at, ngroup = ngroup,
                   seed = seed, trace = trace)

             }

      )

    }

    names(problist) = model.type
    class(problist) = c('hdnom.compare.calibrate')

    problist

  }

#' Print Model Comparison by Calibration Results
#'
#' Print Model Comparison by Calibration Results
#'
#' @param x an object returned by \code{\link{hdnom.compare.calibrate}}.
#' @param ... other parameters (not used).
#'
#' @method print hdnom.compare.calibrate
#'
#' @export
#'
#' @examples
#' NULL
print.hdnom.compare.calibrate = function(x, ...) {

  if (!('hdnom.compare.calibrate' %in% class(x)))
    stop('object class must be "hdnom.compare.calibrate"')

  for (i in 1L:length(x)) {
    print(x[[i]])
    cat('\n')
  }

}

#' Summary of Model Comparison by Calibration Results
#'
#' Summary of Model Comparison by Calibration Results
#'
#' @param object an object returned by \code{\link{hdnom.compare.calibrate}}.
#' @param ... other parameters (not used).
#'
#' @method summary hdnom.compare.calibrate
#'
#' @export
#'
#' @examples
#' NULL
summary.hdnom.compare.calibrate = function(object, ...) {

  if (!('hdnom.compare.calibrate' %in% class(object)))
    stop('object class must be "hdnom.compare.calibrate"')

  for (i in 1L:length(object)) {
    cat('  Model type:', names(object)[i], '\n')
    summary(object[[i]])
    cat('\n')
  }

}

#' Plot Model Comparison by Calibration Results
#'
#' Plot Model Comparison by Calibration Results
#'
#' @param x an object returned by \code{\link{hdnom.compare.calibrate}}.
#' @param xlim x axis limits of the plot.
#' @param ylim y axis limits of the plot.
#' @param col.pal Color palette to use. Possible values are
#' \code{"JCO"}, \code{"Lancet"}, \code{"NPG"}, and \code{"AAAS"}.
#' Default is \code{"JCO"}.
#' @param ... other parameters (not used).
#'
#' @method plot hdnom.compare.calibrate
#'
#' @export
#'
#' @importFrom ggplot2 ggplot aes_string geom_errorbar
#' geom_line geom_point geom_abline scale_colour_manual
#' xlab ylab theme_bw
#'
#' @examples
#' NULL
plot.hdnom.compare.calibrate =
  function(x, xlim = c(0, 1), ylim = c(0, 1),
           col.pal = c('JCO', 'Lancet', 'NPG', 'AAAS'), ...) {

    if (!('hdnom.compare.calibrate' %in% class(x)))
      stop('object class must be "hdnom.compare.calibrate"')

    n = length(x)
    dflist = vector('list', n)

    for (i in 1L:n) {
      dflist[[i]] =
        data.frame('pre' = x[[i]][, 'Predicted'],
                   'obs' = x[[i]][, 'Observed'],
                   'll'  = x[[i]][, 'Lower 95%'],
                   'ul'  = x[[i]][, 'Upper 95%'])
      dflist[[i]][, 'Model'] = names(x)[i]
    }

    df = Reduce('rbind', dflist)

    col.pal = match.arg(col.pal)
    col_pal = switch (
      col.pal,
      JCO   = palette.jco(), Lancet = palette.lancet(),
      NPG   = palette.npg(), AAAS   = palette.aaas())

    ggplot(df, aes_string(x = 'pre', y = 'obs',
                          xmin = xlim[1L], xmax = xlim[2L],
                          ymin = ylim[1L], ymax = ylim[2L],
                          colour = 'Model')) +
      geom_errorbar(aes_string(ymin = 'll', ymax = 'ul'), width = 0.02) +
      geom_line() +
      geom_point(size = 2) +
      geom_abline(slope = 1, intercept = 0, colour = 'grey') +
      xlab('Predicted Survival Probability') +
      ylab('Observed Survival Probability') +
      scale_colour_manual(values = col_pal) +
      theme_bw()

  }
