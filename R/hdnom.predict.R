#' Make Predictions from High-Dimensional Cox Models
#'
#' Predict overall survival probability at certain time points
#' from established Cox models.
#'
#' @param object Model object fitted by \code{hdcox.*()} functions.
#' @param x Data matrix used to fit the model.
#' @param y Response matrix made with \code{\link[survival]{Surv}}.
#' @param newx Matrix (with named columns) of new values for \code{x}
#' at which predictions are to be made.
#' @param pred.at Time point at which prediction should take place.
#' @param ... Other parameters (not used).
#'
#' @return A \code{nrow(newx) x length(pred.at)} matrix containing overall
#' survival probablity.
#'
#' @method predict hdcox.model
#'
#' @importFrom penalized survival
#' @importMethodsFrom penalized predict
#'
#' @export
#'
#' @examples
#' library("survival")
#'
#' # Load imputed SMART data
#' data("smart")
#' x = as.matrix(smart[, -c(1, 2)])
#' time = smart$TEVENT
#' event = smart$EVENT
#' y = Surv(time, event)
#'
#' lassofit = hdcox.lasso(x, y, nfolds = 5, rule = "lambda.1se", seed = 11)
#' predict(lassofit, x, y, newx = x[101:105, ], pred.at = 1:10 * 365)
predict.hdcox.model = function(object, x, y, newx, pred.at, ...) {

  if (!('hdcox.model' %in% class(object)))
    stop('object must be of class "hdcox.model" fitted by hdcox.* functions')

  model.type = gsub('hdcox.model.', '', setdiff(class(object), 'hdcox.model'))

  if (!('matrix' %in% class(newx))) stop('newx must be a matrix')

  time = y[, 1L]
  event = y[, 2L]

  obj.type = switch(model.type,
                    lasso  = 'glmnet', alasso = 'glmnet',
                    enet   = 'glmnet', aenet  = 'glmnet',
                    mcp    = 'ncvreg', mnet   = 'ncvreg',
                    scad   = 'ncvreg', snet   = 'ncvreg',
                    flasso = 'penalized'
  )

  switch(obj.type,
         glmnet = {

           lp = predict(object[[paste0(model.type, '_model')]], x, type = 'link')
           basesurv = glmnet.basesurv(time, event, lp, pred.at)
           lpnew = predict(object[[paste0(model.type, '_model')]], newx, type = 'link')
           p = exp(exp(lpnew) %*% -t(basesurv$'cumulative_base_hazard'))

         },
         ncvreg = {

           lp = predict(object[[paste0(model.type, '_model')]], x, type = 'link')
           basesurv = ncvreg.basesurv(time, event, lp, pred.at)
           lpnew = predict(object[[paste0(model.type, '_model')]], newx, type = 'link')
           p = exp(exp(lpnew) %*% -t(basesurv$'cumulative_base_hazard'))

           # # alternative method using ncvreg built-in prediction directly
           # # almost identical results, but sometimes produces NAs in prediction
           # # e.g. pred.at = 1:10 * 365
           # survfun = predict(object[[paste0(model.type, '_model')]], newx, type = 'survival')
           # p = matrix(NA, nrow = nrow(newx), ncol = length(pred.at))
           # for (i in 1L:nrow(newx)) p[i, ] = sapply(pred.at, survfun[[i]])

         },
         penalized = {

           pred = predict(object[[paste0(model.type, '_model')]], newx)
           p = matrix(NA, nrow = nrow(newx), ncol = length(pred.at))
           for (i in 1L:length(pred.at)) p[, i] = survival(pred, time = pred.at[i])

         }
  )

  colnames(p) = as.character(pred.at)
  return(p)

}

#' Extract Information of Selected Variables from High-Dimensional Cox Models
#'
#' Extract the names and type of selected variables from established
#' high-dimensional Cox models.
#'
#' @param object Model object fitted by \code{hdcox.*()} functions.
#' @param x Data matrix used to fit the model.
#'
#' @export selected.info
#'
#' @examples
#' NULL
selected.info = function(object, x) {

  # type (logical or continous)
  #
  # if logical:
  # return:
  # two possible values
  # (select box)
  #
  # if categorical (more than 2 values, all int)
  # return:
  # lower bound
  # upper bound
  # (numericInput)
  #
  # if continuous (non-int detected):
  # return:
  # upper bound
  # lower bound
  # (textinput)

  varinfo = NULL
  return(varinfo)

}
