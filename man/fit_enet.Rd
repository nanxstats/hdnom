% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/1_1_model.R
\name{fit_enet}
\alias{fit_enet}
\title{Model selection for high-dimensional Cox models with elastic-net penalty}
\usage{
fit_enet(
  x,
  y,
  nfolds = 5L,
  alphas = seq(0.05, 0.95, 0.05),
  rule = c("lambda.min", "lambda.1se"),
  seed = 1001,
  parallel = FALSE
)
}
\arguments{
\item{x}{Data matrix.}

\item{y}{Response matrix made by \code{\link[survival]{Surv}}.}

\item{nfolds}{Fold numbers of cross-validation.}

\item{alphas}{Alphas to tune in \code{\link[glmnet]{cv.glmnet}}.}

\item{rule}{Model selection criterion, \code{"lambda.min"} or
\code{"lambda.1se"}. See \code{\link[glmnet]{cv.glmnet}}
for details.}

\item{seed}{A random seed for cross-validation fold division.}

\item{parallel}{Logical. Enable parallel parameter tuning or not,
default is \code{FALSE}. To enable parallel tuning, load the
\code{doParallel} package and run \code{registerDoParallel()}
with the number of CPU cores before calling this function.}
}
\description{
Automatic model selection for high-dimensional Cox models
with elastic-net penalty, evaluated by penalized partial-likelihood.
}
\examples{
data("smart")
x <- as.matrix(smart[, -c(1, 2)])
time <- smart$TEVENT
event <- smart$EVENT
y <- survival::Surv(time, event)

# To enable parallel parameter tuning, first run:
# library("doParallel")
# registerDoParallel(detectCores())
# then set fit_enet(..., parallel = TRUE).

fit <- fit_enet(
  x, y,
  nfolds = 3, alphas = c(0.3, 0.7),
  rule = "lambda.min", seed = 3
)

nom <- as_nomogram(
  fit, x, time, event,
  pred.at = 365 * 2,
  funlabel = "2-Year Overall Survival Probability"
)

plot(nom)
}
