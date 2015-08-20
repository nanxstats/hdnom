#' Nomograms for Cox models fitted with glmnet
#'
#' Nomograms for Cox models fitted with glmnet
#'
#' @param Fitted \code{"glmnet"} model object.
#' @param x Matrix of training data used for the \code{glmnet} object.
#' @param ddist Data frame version of x, made by \code{datadist()}
#' @param s Value of the penalty parameter lambda in \code{glmnet}.
#' We will use the selected variables at the provided \code{s} to
#' build the nomogram, and make predictions.
#' See the example for choosing a proper lambda value extracted
#' from cross-validation results.
#' @param ... other arguments for \code{nomogram}.
#'
#' @export glmnet.nomogram
#'
#' @importFrom rms ols nomogram
#'
#' @examples
#' TBA
glmnet.nomogram = function(object, x, ddist, s, pred.at, funlabel) {

  if (!all(c('coxnet', 'glmnet') %in% class(object)))
    stop('object class must be "glmnet" and "coxnet"')

  glmnet_pred_lp = as.vector(predict(object, newx = x, s = s, type = 'link'))

  all_vars = rownames(object$beta)
  selected_vars = all_vars[which(abs(coef(fit, s = s)) > .Machine$double.eps)]
  ols_formula = paste('glmnet_pred_lp ~',
                      paste(selected_vars, collapse = ' + '))
  ols_fit = ols(as.formula(ols_formula), data = ddist,
                sigma = 1, x = TRUE, y = TRUE)

  nom = list('ols_fit' = ols_fit, 'pred.at' = pred.at, 'funlabel' = funlabel)

  class(nom) = 'glmnet.nomogram'

  nom

}

plot.glmnet.nomogram = function (object) {

  if (class(object) != 'glmnet.nomogram')
    stop('object class must be "glmnet.nomogram"')

  nom = nomogram(fit = object$ols_fit)
  # TODO: Arguments cannot be passed to here ...
  # investigate svycox.nomogram to do this only allowing Overall Survival Prob.

  plot(nom)

}

library('glmnet')
library('survival')
library('rms')

# Original SMART data not truncated
smart = read.table('~/Desktop/hdnomo-support/smartdata/SMARTs.tsv', header = TRUE, sep = '\t')
smart = na.omit(smart)
x = as.matrix(smart[, -c(1, 2)])
time = smart[, 1]
status = smart[, 2]

x.df = as.data.frame(x)
dd = datadist(x.df)
options(datadist = "dd")

set.seed(1010)
cvfit = cv.glmnet(x, Surv(time, status), family = "cox", nfolds = 5)
fit = glmnet(x, Surv(time, status), family = "cox")
nom = glmnet.nomogram(fit, x, ddist = x.df, s = cvfit$lambda.min,
                      pred.at = 365, funlabel = '2-Year Overall Survival Probability')
plot(nom)



library('peperr')
?predictProb

library('c060')
predictProb.glmnet()

library('SvyNom')
?svycox.nomogram
library('survey')
library('rms')
data(noNA)
dd=datadist(noNA)
options(datadist="dd")
dstr2=svydesign(id=~1, strata=~group, prob=~inv_weight, fpc=~ssize, data=noNA)
mynom=svycox.nomogram(.design=dstr2,
                      .model=Surv(survival,surv_cens)~ECOG+liver_only+Alb+Hb+Age+Differentiation + Gt_1_m1site+lymph_only,
                      .data=noNA, pred.at=24, fun.lab="Prob of 2 Yr OS")
plot(mynom$nomog)
.design=dstr2
.model=Surv(survival,surv_cens)~ECOG+liver_only+Alb+Hb+Age+Differentiation + Gt_1_m1site+lymph_only
.data=noNA
pred.at=24
fun.lab="Prob of 2 Yr OS"



# this is ok
time.at = pred.survey.cox[[1]]$time[which(pred.survey.cox[[1]]$time > pred.at)[1] - 1]

length(pred.survey.cox)
pred.survey.cox[[10]]$time == pred.survey.cox[[20]]$time
slength(pred.survey.cox[[20]]$time)
oneidx = which(noNA$surv_cens==1)
length(oneidx)
survtime = noNA[oneidx, 'survival']
names(survtime) = oneidx
all(sort(survtime) == pred.survey.cox[[10]]$time)

# this is TBD
.baseline = exp(log(pred.survey.cox[[1]]$surv[names(time.at)])/exp(svy.cox.fit$linear.predictors[1]))

str(pred.survey.cox[[1]])

oneidx.glmnet = which(status == 1)
survtime.glmnet = time[oneidx.glmnet]
names(survtime.glmnet) = oneidx.glmnet
survtime.glmnet = sort(survtime.glmnet)

library('c060')
predictProb(object = fit, response = Surv(time, status),
            x = x, times = survtime.glmnet,
            complexity = cvfit$lambda.min)  # oneidx.glmnet should be sorted!

z = Surv(time, status)
?predictProb.coxnet

range(pred.survey.cox[[1]]$surv)
svy.cox.fit$linear.predictors

.tempfun = function(x) .baseline[[1]]^exp(x)





## dig into c060::predictProb.coxnet for estimating baseline hazard
fit.glmnet = function (response, x, cplx, ...) {

  res = NULL
  tryerr = try(res <- glmnet(y = response, x = data.matrix(x), lambda = cplx,  ...), silent = TRUE)

  if (class(tryerr) != 'try-error' && 'coxnet' %in% class(res)) {
    res$linear.predictor = as.numeric(predict(res, newx=data.matrix(x), type="link"))
    res$response = response
  }

  res

}

predictProb.coxnet = predictProb.glmnet = function (object, response, x,
                                                    times, complexity,  ...) {

  lp = as.numeric(predict(object, newx = data.matrix(x), s = complexity, type = 'link'))
  basesurv = basesurv(object$response, object$linear.predictor, sort(unique(times)))
  p = exp(exp(lp) %*% -t(basesurv$cumBaseHaz))

  if (NROW(p) != NROW(x) || NCOL(p) != length(times))
    stop("Prediction failed")

  p

}

# baseline survival/hazard Breslow estimator
# function essentially based on gbm::basehaz.gbm
basesurv = function (response, lp, times.eval = NULL, centered = FALSE) {

  if (is.null(times.eval)) times.eval = sort(unique(response[, 1]))

  t.unique = sort(unique(response[, 1][response[, 2] == 1]))
  alpha = length(t.unique)

  for (i in 1:length(t.unique)) {
    alpha[i] = sum(response[, 1][response[, 2] == 1] ==
                     t.unique[i])/sum(exp(lp[response[, 1] >=  t.unique[i]]))
  }

  obj = approx(t.unique, cumsum(alpha), yleft = 0, xout = times.eval, rule = 2)

  if (centered) obj$y = obj$y * exp(mean(lp))
  obj$z = exp(-obj$y)

  names(obj) = c('times', 'cumBaseHaz', 'BaseSurv')

  obj

}

library('gbm')
?basehaz.gbm

basehaz.gbm =  function (t, delta, f.x, t.eval = NULL, smooth = FALSE, cumulative = TRUE) {
  t.unique = sort(unique(t[delta == 1]))
  alpha = length(t.unique)
  for (i in 1:length(t.unique)) {
    alpha[i] = sum(t[delta == 1] == t.unique[i])/sum(exp(f.x[t >= t.unique[i]]))
  }
  if (!smooth && !cumulative) {
    if (!is.null(t.eval)) {
      stop("Cannot evaluate unsmoothed baseline hazard at t.eval.")
    }
  }
  else if (smooth && !cumulative) {
    lambda.smooth = supsmu(t.unique, alpha)
  }
  else if (smooth && cumulative) {
    lambda.smooth = supsmu(t.unique, cumsum(alpha))
  }
  else {
    lambda.smooth = list(x = t.unique, y = cumsum(alpha))
  }
  if (!is.null(t.eval)) {
    obj = approx(lambda.smooth$x, lambda.smooth$y, xout = t.eval)$y
  }
  else {
    obj = approx(lambda.smooth$x, lambda.smooth$y, xout = t)$y
  }

  obj

}
