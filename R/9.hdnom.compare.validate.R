#' Compare High-Dimensional Cox Models by Model Validation
#'
#' Compare High-Dimensional Cox Models by Model Validation
#'
#' @param x Matrix of training data used for fitting the model;
#' on which to run the validation.
#' @param time Survival time.
#' Must be of the same length with the number of rows as \code{x}.
#' @param event Status indicator, normally 0 = alive, 1 = dead.
#' Must be of the same length with the number of rows as \code{x}.
#' @param model.type Model types to compare. Could be at least two
#' of \code{"lasso"}, \code{"alasso"}, \code{"flasso"}, \code{"enet"},
#' \code{"aenet"}, \code{"mcp"}, \code{"mnet"}, \code{"scad"},
#' or \code{"snet"}.
#' @param method Validation method.
#' Could be \code{"bootstrap"}, \code{"cv"}, or \code{"repeated.cv"}.
#' @param boot.times Number of repetitions for bootstrap.
#' @param nfolds Number of folds for cross-validation and
#' repeated cross-validation.
#' @param rep.times Number of repeated times for repeated cross-validation.
#' @param tauc.type Type of time-dependent AUC.
#' Including \code{"CD"} proposed by Chambless and Diao (2006).,
#' \code{"SZ"} proposed by Song and Zhou (2008).,
#' \code{"UNO"} proposed by Uno et al. (2007).
#' @param tauc.time Numeric vector. Time points at which to evaluate
#' the time-dependent AUC.
#' @param seed A random seed for cross-validation fold division.
#' @param trace Logical. Output the validation progress or not.
#' Default is \code{TRUE}.
#'
#' @export hdnom.compare.validate
#'
#' @importFrom ggplot2 ggplot
#'
#' @references
#' Chambless, L. E. and G. Diao (2006).
#' Estimation of time-dependent area under the ROC curve for long-term
#' risk prediction.
#' \emph{Statistics in Medicine} 25, 3474--3486.
#'
#' Song, X. and X.-H. Zhou (2008).
#' A semiparametric approach for the covariate specific ROC curve with
#' survival outcome.
#' \emph{Statistica Sinica} 18, 947--965.
#'
#' Uno, H., T. Cai, L. Tian, and L. J. Wei (2007).
#' Evaluating prediction rules for t-year survivors with censored
#' regression models.
#' \emph{Journal of the American Statistical Association} 102, 527--537.
#'
#' @examples
#' # Load imputed SMART data
#' data(smart)
#' x = as.matrix(smart[, -c(1, 2)])[1:1000, ]
#' time = smart$TEVENT[1:1000]
#' event = smart$EVENT[1:1000]
#'
#' # Compare lasso and adaptive lasso by 5-fold cross-validation
#' cmp.val.cv =
#'   hdnom.compare.validate(x, time, event,
#'                          model.type = c("lasso", "alasso"),
#'                          method = "cv", nfolds = 5, tauc.type = "UNO",
#'                          tauc.time = seq(0.25, 2, 0.25) * 365, seed = 1001)
#'
#' print(cmp.val.cv)
#' summary(cmp.val.cv)
#' plot(cmp.val.cv)
#' plot(cmp.val.cv, interval = TRUE)
hdnom.compare.validate =
  function(x, time, event,
           model.type = c('lasso', 'alasso', 'flasso',
                          'enet', 'aenet',
                          'mcp', 'mnet',
                          'scad', 'snet'),
           method = c('bootstrap', 'cv', 'repeated.cv'),
           boot.times = NULL, nfolds = NULL, rep.times = NULL,
           tauc.type = c("CD", "SZ", "UNO"), tauc.time,
           seed = 1001, trace = TRUE) {

    method = match.arg(method)
    tauc.type = match.arg(tauc.type)

    if (!all(model.type %in% c('lasso', 'alasso', 'flasso', 'enet', 'aenet',
                               'mcp', 'mnet', 'scad', 'snet'))) {
      stop('Unknown model type(s) specified')
    }

    nmodel = length(model.type)
    tauclist = vector('list', nmodel)

    # check parameters for different methods
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

               tauclist[[i]] =
                 hdnom.validate(
                   x, time, event, model.type = 'lasso',
                   alpha = 1, lambda = cvfit$'lasso_best_lambda',
                   method = method, boot.times = boot.times, nfolds = nfolds, rep.times = rep.times,
                   tauc.type = tauc.type, tauc.time = tauc.time,
                   seed = seed, trace = trace)

             },

             alasso = {

               cvfit = hdcox.alasso(x, Surv(time, event), nfolds = 5L,
                                    rule = 'lambda.1se', seed = rep(seed, 2))

               tauclist[[i]] =
                 hdnom.validate(
                   x, time, event, model.type = 'alasso',
                   alpha = 1, lambda = cvfit$'alasso_best_lambda', pen.factor = cvfit$'pen_factor',
                   method = method, boot.times = boot.times, nfolds = nfolds, rep.times = rep.times,
                   tauc.type = tauc.type, tauc.time = tauc.time,
                   seed = seed, trace = trace)

             },

             flasso = {

               cvfit = hdcox.flasso(x, Surv(time, event), nfolds = 5L, seed = seed)

               tauclist[[i]] =
                 hdnom.validate(
                   x, time, event, model.type = 'flasso',
                   alpha = 1, lambda = cvfit$'flasso_best_lambda',
                   method = method, boot.times = boot.times, nfolds = nfolds, rep.times = rep.times,
                   tauc.type = tauc.type, tauc.time = tauc.time,
                   seed = seed, trace = trace)

             },

             enet = {

               cvfit = hdcox.enet(
                 x, Surv(time, event), nfolds = 5L,
                 alphas = c(0.1, 0.25, 0.5, 0.75, 0.9),  # to reduce computation time
                 rule = 'lambda.1se', seed = seed)

               tauclist[[i]] =
                 hdnom.validate(
                   x, time, event, model.type = 'enet',
                   alpha = cvfit$'enet_best_alpha', lambda = cvfit$'enet_best_lambda',
                   method = method, boot.times = boot.times, nfolds = nfolds, rep.times = rep.times,
                   tauc.type = tauc.type, tauc.time = tauc.time,
                   seed = seed, trace = trace)

             },

             aenet = {

               cvfit = hdcox.aenet(
                 x, Surv(time, event), nfolds = 5L,
                 alphas = c(0.1, 0.25, 0.5, 0.75, 0.9),  # to reduce computation time
                 rule = 'lambda.1se', seed = rep(seed, 2))

               tauclist[[i]] =
                 hdnom.validate(
                   x, time, event, model.type = 'aenet',
                   alpha = cvfit$'aenet_best_alpha', lambda = cvfit$'aenet_best_lambda', pen.factor = cvfit$'pen_factor',
                   method = method, boot.times = boot.times, nfolds = nfolds, rep.times = rep.times,
                   tauc.type = tauc.type, tauc.time = tauc.time,
                   seed = seed, trace = trace)

             },

             mcp = {

               cvfit = hdcox.mcp(x, Surv(time, event), nfolds = 5L, seed = seed)

               tauclist[[i]] =
                 hdnom.validate(
                   x, time, event, model.type = 'mcp',
                   alpha = 1, gamma = cvfit$'mcp_best_gamma', lambda = cvfit$'mcp_best_lambda',
                   method = method, boot.times = boot.times, nfolds = nfolds, rep.times = rep.times,
                   tauc.type = tauc.type, tauc.time = tauc.time,
                   seed = seed, trace = trace)

             },

             mnet = {

               cvfit = hdcox.mnet(
                 x, Surv(time, event), nfolds = 5L,
                 alphas = c(0.1, 0.25, 0.5, 0.75, 0.9), # to reduce computation time
                 seed = seed)

               tauclist[[i]] =
                 hdnom.validate(
                   x, time, event, model.type = 'mnet',
                   alpha = cvfit$'mnet_best_alpha', gamma = cvfit$'mnet_best_gamma', lambda = cvfit$'mnet_best_lambda',
                   method = method, boot.times = boot.times, nfolds = nfolds, rep.times = rep.times,
                   tauc.type = tauc.type, tauc.time = tauc.time,
                   seed = seed, trace = trace)

             },

             scad = {

               cvfit = hdcox.scad(x, Surv(time, event), nfolds = 5L, seed = seed)

               tauclist[[i]] =
                 hdnom.validate(
                   x, time, event, model.type = 'scad',
                   alpha = 1, gamma = cvfit$'scad_best_gamma', lambda = cvfit$'scad_best_lambda',
                   method = method, boot.times = boot.times, nfolds = nfolds, rep.times = rep.times,
                   tauc.type = tauc.type, tauc.time = tauc.time,
                   seed = seed, trace = trace)

             },

             snet = {

               cvfit = hdcox.snet(
                 x, Surv(time, event), nfolds = 5L,
                 alphas = c(0.1, 0.25, 0.5, 0.75, 0.9), # to reduce computation time
                 seed = seed)

               tauclist[[i]] =
                 hdnom.validate(
                   x, time, event, model.type = 'snet',
                   alpha = cvfit$'snet_best_alpha', gamma = cvfit$'snet_best_gamma', lambda = cvfit$'snet_best_lambda',
                   method = method, boot.times = boot.times, nfolds = nfolds, rep.times = rep.times,
                   tauc.type = tauc.type, tauc.time = tauc.time,
                   seed = seed, trace = trace)

             }

      )

    }

    names(tauclist) = model.type
    class(tauclist) = c('hdnom.compare.validate')

    tauclist

  }

#' Print Model Comparison by Validation Results
#'
#' Print Model Comparison by Validation Results
#'
#' @param x An object returned by \code{\link{hdnom.compare.validate}}.
#' @param ... Other parameters (not used).
#'
#' @method print hdnom.compare.validate
#'
#' @export
#'
#' @examples
#' NULL
print.hdnom.compare.validate = function(x, ...) {

  if (!('hdnom.compare.validate' %in% class(x)))
    stop('object class must be "hdnom.compare.validate"')

  for (i in 1L:length(x)) {
    print(x[[i]])
    cat('\n\n')
  }

}

#' Summary of Model Comparison by Validation Results
#'
#' Summary of Model Comparison by Validation Results
#'
#' @param object An object \code{\link{hdnom.compare.validate}}.
#' @param silent Print summary table header or not,
#' default is \code{FALSE}.
#' @param ... Other parameters (not used).
#'
#' @method summary hdnom.compare.validate
#'
#' @export
#'
#' @examples
#' NULL
summary.hdnom.compare.validate = function(object, silent = FALSE, ...) {

  if (!('hdnom.compare.validate' %in% class(object)))
    stop('object class must be "hdnom.compare.validate"')

  for (i in 1L:length(object)) {
    cat('Model type:', names(object)[i], '\n')
    print(summary(object[[i]], silent = TRUE))
    cat('\n')
  }

}

#' Plot Model Comparison by Validation Results
#'
#' Plot Model Comparison by Validation Results
#'
#' @param x An object returned by \code{\link{hdnom.compare.validate}}.
#' @param interval Show maximum, minimum, 0.25, and 0.75 quantiles of
#' time-dependent AUC as ribbons? Default is \code{FALSE}.
#' @param col.pal Color palette to use. Possible values are
#' \code{"JCO"}, \code{"Lancet"}, \code{"NPG"}, and \code{"AAAS"}.
#' Default is \code{"JCO"}.
#' @param ... Other parameters (not used).
#'
#' @method plot hdnom.compare.validate
#'
#' @importFrom ggplot2 ggplot aes_string geom_point geom_line
#' scale_x_continuous scale_colour_manual theme_bw ylab
#'
#' @export
#'
#' @examples
#' NULL
plot.hdnom.compare.validate =
  function(x, interval = FALSE,
           col.pal = c('JCO', 'Lancet', 'NPG', 'AAAS'), ...) {

    if (!('hdnom.compare.validate' %in% class(x)))
      stop('object class must be "hdnom.compare.validate"')

    n = length(x)
    dflist = vector('list', n)

    for (i in 1L:n) {
      dflist[[i]] = as.data.frame(t(summary(x[[i]], silent = TRUE)))
      tauc_time = attr(x[[i]], 'tauc.time')

      # special processing for repeated cv
      if (any(grepl(pattern = 'validate.repeated.cv', class(x[[i]]))))
        names(dflist[[i]]) = sapply(strsplit(names(dflist[[i]]), 'Mean of '), '[', 2L)

      dflist[[i]][, 'Time'] = tauc_time
      dflist[[i]][, 'Model'] = names(x)[i]
      names(dflist[[i]])[which(names(dflist[[i]]) == '0.25 Qt.')] = 'Qt25'
      names(dflist[[i]])[which(names(dflist[[i]]) == '0.75 Qt.')] = 'Qt75'
    }

    df = Reduce('rbind', dflist)

    col.pal = match.arg(col.pal)
    col_pal = switch (
      col.pal,
      JCO   = palette.jco(), Lancet = palette.lancet(),
      NPG   = palette.npg(), AAAS   = palette.aaas())

    if (!interval) {

      ggplot(data = df, aes_string(x = 'Time', y = 'Mean',
                                   colour = 'Model', fill = 'Model')) +
        geom_point() +
        geom_line() +
        geom_point(data = df, aes_string(x = 'Time', y = 'Median',
                                         colour = 'Model', fill = 'Model')) +
        geom_line(data = df, aes_string(x = 'Time', y = 'Median',
                                        colour = 'Model', fill = 'Model'),
                  linetype = 'dashed') +
        scale_x_continuous(breaks = df$'Time') +
        scale_colour_manual(values = col_pal) +
        theme_bw() +
        ylab('Area under ROC')

    } else {

      ggplot(data = df, aes_string(x = 'Time', y = 'Mean',
                                   colour = 'Model', fill = 'Model')) +
        geom_point() +
        geom_line() +
        geom_point(data = df, aes_string(x = 'Time', y = 'Median',
                                         colour = 'Model', fill = 'Model')) +
        geom_line(data = df, aes_string(x = 'Time', y = 'Median',
                                        colour = 'Model', fill = 'Model'),
                  linetype = 'dashed') +
        geom_ribbon(data = df, aes_string(ymin = 'Qt25', ymax = 'Qt75',
                                          colour = 'Model', fill = 'Model'),
                    linetype = 0, alpha = 0.1) +
        geom_ribbon(data = df, aes_string(ymin = 'Min', ymax = 'Max',
                                          colour = 'Model', fill = 'Model'),
                    linetype = 0, alpha = 0.05) +
        scale_x_continuous(breaks = df$'Time') +
        scale_colour_manual(values = col_pal) +
        scale_fill_manual(values = col_pal) +
        theme_bw() +
        ylab('Area under ROC')

    }

  }
