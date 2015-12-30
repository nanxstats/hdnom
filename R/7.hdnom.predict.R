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
           # # almost identical results, but sometimes produces NAs in practice
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
#' @export hdnom.varinfo
#'
#' @return A list containing the name, type and range of the selected variables.
#'
#' @examples
#' library("glmnet")
#' library("survival")
#'
#' # Load imputed SMART data
#' data("smart")
#' x = as.matrix(smart[, -c(1, 2)])
#' time = smart$TEVENT
#' event = smart$EVENT
#' y = Surv(time, event)
#'
#' # Fit Cox model by lasso penalization
#' lassofit = hdcox.lasso(x, y, nfolds = 5, rule = "lambda.1se", seed = 11)
#' hdnom.varinfo(lassofit, x)
hdnom.varinfo = function(object, x) {

  if (!('hdcox.model' %in% class(object)))
    stop('object must be of class "hdcox.model" fitted by hdcox.* functions')

  model.type = gsub('hdcox.model.', '', setdiff(class(object), 'hdcox.model'))

  obj.type = switch(model.type,
                    lasso  = 'glmnet', alasso = 'glmnet',
                    enet   = 'glmnet', aenet  = 'glmnet',
                    mcp    = 'ncvreg', mnet   = 'ncvreg',
                    scad   = 'ncvreg', snet   = 'ncvreg',
                    flasso = 'penalized'
  )

  switch(obj.type,
         glmnet = {
           nonzero_idx = which(as.logical(abs(object[[paste0(model.type, '_model')]][['beta']]) > .Machine$double.eps))
           nonzero_var = rownames(object[[paste0(model.type, '_model')]][['beta']])[nonzero_idx]
         },
         ncvreg = {
           nonzero_idx = which(object[[paste0(model.type, '_model')]][['beta']][-1L, ] > .Machine$double.eps)
           nonzero_var = names(object[[paste0(model.type, '_model')]][['beta']][-1L, ])[nonzero_idx]
         },
         penalized = {
           nonzero_idx = which(object[[paste0(model.type, '_model')]]@'penalized' > .Machine$double.eps)
           nonzero_var = colnames(x)[nonzero_idx]
         }
  )

  varinfo = list('name' = NULL, 'type' = NULL, 'domain' = NULL)

  varinfo[['name']] = nonzero_var
  nvar = length(varinfo[['name']])

  is.wholenumber = function(x, tol = .Machine$double.eps^0.5) abs(x - round(x)) < tol

  # type: logical, categorical or continuous
  for (i in 1L:nvar) {
    var = x[, varinfo[['name']][i]]
    if (all(is.wholenumber(var)) & nlevels(as.factor(var)) == 2L) {
      varinfo[['type']][i] = 'logical'
      varinfo[['domain']][[i]] = unique(var)
    } else if (all(is.wholenumber(var)) & nlevels(as.factor(var)) > 2L) {
      varinfo[['type']][i] = 'categorical'
      varinfo[['domain']][[i]] = c(min(var), max(var))
    } else if (any(!is.wholenumber(var))) {
      varinfo[['type']][i] = 'continuous'
      varinfo[['domain']][[i]] = c(min(var), max(var))
    } else {
      stop(paste0('unrecognized variable type: ', varinfo[['name']][i]))
    }
  }

  class(varinfo) = 'hdnom.varinfo'
  return(varinfo)

}
