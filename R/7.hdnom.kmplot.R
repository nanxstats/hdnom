#' Kaplan-Meier Plot with Number at Risk Table for Internal Calibration and
#' External Calibration Results
#'
#' Kaplan-Meier Plot with Number at Risk Table for Internal Calibration and
#' External Calibration Results
#'
#' @param object An object returned by \code{\link{hdnom.calibrate}} or
#' \code{\link{hdnom.external.calibrate}}.
#' @param group.name Risk group labels. Default is
#' Group 1, Group 2, ..., Group k.
#' @param time.at Time points to evaluate the number at risk.
#' @param col.pal Color palette to use. Possible values are
#' \code{"JCO"}, \code{"Lancet"}, \code{"NPG"}, and \code{"AAAS"}.
#' Default is \code{"JCO"}.
#'
#' @importFrom survival survfit
#' @importFrom survival Surv
#' @importFrom stats formula
#'
#' @export hdnom.kmplot
#'
#' @examples
#' library("survival")
#'
#' # Load imputed SMART data
#' data("smart")
#'
#' # Use the first 1000 samples as training data
#' # (the data used for internal validation)
#' x = as.matrix(smart[, -c(1, 2)])[1:1000, ]
#' time = smart$TEVENT[1:1000]
#' event = smart$EVENT[1:1000]
#'
#' # Take the next 1000 samples as external calibration data
#' # In practice, usually use data collected in other studies
#' x_new = as.matrix(smart[, -c(1, 2)])[1001:2000, ]
#' time_new = smart$TEVENT[1001:2000]
#' event_new = smart$EVENT[1001:2000]
#'
#' # Fit Cox model by lasso penalization
#' lassofit = hdcox.lasso(x, Surv(time, event), nfolds = 5, rule = "lambda.1se", seed = 11)
#'
#' ### Internal calibration
#' cal.int = hdnom.calibrate(x, time, event, model.type = "lasso",
#'                           alpha = 1, lambda = lassofit$'lasso_best_lambda',
#'                           method = "cv", nfolds = 5,
#'                           pred.at = 365 * 9, ngroup = 3)
#'
#' hdnom.kmplot(cal.int, group.name = c('High risk', 'Medium risk', 'Low risk'),
#'              time.at = 1:6 * 365)
#'
#' ### External calibration
#' cal.ext =
#'   hdnom.external.calibrate(lassofit, x, time, event,
#'                            x_new, time_new, event_new,
#'                            pred.at = 365 * 5, ngroup = 3)
#'
#' hdnom.kmplot(cal.ext, group.name = c('High risk','Medium risk', 'Low risk'),
#'              time.at = 1:6 * 365)
hdnom.kmplot = function(object, group.name = NULL, time.at = NULL,
                        col.pal = c('JCO', 'Lancet', 'NPG', 'AAAS')) {

  if (!(any(c('hdnom.calibrate', 'hdnom.external.calibrate') %in% class(object))))
    stop('object class must be "hdnom.calibrate" or "hdnom.external.calibrate"')

  time = attr(object, 'surv.time')
  event = attr(object, 'surv.event')
  grp = attr(object, 'risk.group')
  df = data.frame(time, event, grp)
  fit = survfit(formula('Surv(time, event) ~ grp'))
  col.pal = match.arg(col.pal)
  kmplot(fit = fit, group.name = group.name,
         time.at = time.at, surv.df = df, col.pal = col.pal)

}

#' Kaplan-Meier Plot with Number at Risk Table
#'
#' @param fit a \code{\link[survival]{survfit}} object
#' @param group.name Group labels. Default is Group 1, Group 2, ... Group k.
#' @param time.at Time points to evaluate the number at risk.
#' @param surv.df Data frame containing survival time, event and risk group
#' for log-rank test.
#' @param col.pal color palette to use
#'
#' @importFrom stats pchisq
#' @importFrom survival survdiff
#' @importFrom gridExtra grid.arrange arrangeGrob
#' @importFrom ggplot2 geom_step geom_blank geom_text element_blank
#' element_line element_text element_rect scale_colour_manual
#' scale_y_discrete unit annotate
#'
#' @keywords internal
kmplot = function(fit, group.name = NULL, time.at = NULL,
                  surv.df = NULL, col.pal = NULL) {

  if (is.null(group.name))
    group.name = paste('Group', gsub('grp=', '', levels(summary(fit)$'strata')))
  if (is.null(time.at)) stop('time.at must be specified')
  if (is.null(surv.df)) stop('surv.df must be specified')

  # kaplan-meier data
  km_df =
    data.frame(time = fit$'time', surv = fit$'surv',
               upper = fit$'upper', lower = fit$'lower',
               n.risk = fit$'n.risk', n.event = fit$'n.event',
               risk.group = summary(fit, censored = TRUE)$'strata')
  levels(km_df$'risk.group') = group.name
  zero_df =
    data.frame(time = 0, surv = 1,
               upper = 1, lower = 1,
               n.risk = NA, n.event = NA,
               risk.group = factor(group.name, levels = levels(km_df$'risk.group')))
  km_df = rbind(zero_df, km_df)

  if (is.null(col.pal)) stop('col.pal must be specified')

  col_pal = switch (
    col.pal,
    JCO   = palette.jco(), Lancet = palette.lancet(),
    NPG   = palette.npg(), AAAS   = palette.aaas())

  # kaplan-meier plot
  k = max(nchar(group.name))
  km_plot = ggplot(km_df,
                   aes_string(x = 'time', y = 'surv', group = 'risk.group')) +
    geom_step(aes_string(colour = 'risk.group'), size = 0.7) +
    scale_x_continuous(breaks = time.at) +
    scale_colour_manual(values = col_pal) +
    theme_bw() +
    theme(panel.grid.major = element_blank(),  # remove grid lines + half-open frame
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          panel.border = element_blank(),
          # axis.line = element_line(colour = 'black')
          # no axis shown due to a bug introduced in ggplot2 2.1.0
          # workaround: replaced by the two lines below
          axis.line.x = element_line(colour = 'black'),
          axis.line.y = element_line(colour = 'black')) +
    theme(axis.title.x = element_text(vjust = 0.5)) +
    theme(legend.key = element_rect(colour = NA),
          legend.title = element_blank()) +
    theme(plot.margin = unit(c(0, 1, 0.5, ifelse(k < 10, 1.5, 2.5)), 'lines')) +
    xlab('Time') +
    ylab('Overall Survival Probability')

  # add log-rank test p-value
  sdiff = survdiff(formula('Surv(time, event) ~ grp'), data = surv.df)
  pval = pchisq(sdiff$'chisq', length(sdiff$'n') - 1L, lower.tail = FALSE)
  text = ifelse(pval < 0.001, 'Log-rank P < 0.001', paste('Log-rank P =', signif(pval, 3)))
  km_plot = km_plot + annotate('text', x = 0.1 * max(fit$'time'), y = 0, label = text)

  ## place holder plot
  placeholder_plot = ggplot(km_df, aes_string(x = 'time', y = 'surv')) +
    geom_blank() +
    theme_bw() +
    theme(axis.text.x = element_blank(), axis.text.y = element_blank(),
          axis.title.x = element_blank(), axis.title.y = element_blank(),
          axis.ticks = element_blank(), panel.grid.major = element_blank(),
          panel.border = element_blank())

  # number at risk data
  table_df =
    data.frame(time = summary(fit, times = time.at, extend = TRUE)$'time',
               risk.group = summary(fit, times = time.at, extend = TRUE)$'strata',
               n.risk = summary(fit, times = time.at, extend = TRUE)$'n.risk')

  # number at risk plot
  table_plot = ggplot(table_df,
                      aes_string(x = 'time', y = 'risk.group',
                                 label = format('n.risk', nsmall = 0))) +
    geom_text(size = 4) +
    scale_y_discrete(breaks = as.character(levels(table_df$'risk.group')),
                     labels = group.name) +
    scale_x_continuous('Number at risk', limits = c(0, max(fit$'time'))) +
    theme_bw() +
    theme(axis.title.x = element_text(vjust = 1),
          axis.text.y = element_text(hjust = 1),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks = element_blank()) +
    theme(plot.margin = unit(c(-1.5, 1, 0.1, ifelse(k < 10, 2.5, 3.5) - 0.28 * k), 'lines')) +
    xlab(NULL) +
    ylab(NULL)

  grid.arrange(km_plot, placeholder_plot, table_plot,
               nrow = 3, ncol = 1, clip = FALSE,
               heights = unit(c(2, 0.1, 0.25), c('null', 'null', 'null')))

  p = arrangeGrob(km_plot, placeholder_plot, table_plot,
                  nrow = 3, ncol = 1, clip = FALSE,
                  heights = unit(c(2, 0.1, 0.25), c('null', 'null', 'null')))

  invisible(p)

}
