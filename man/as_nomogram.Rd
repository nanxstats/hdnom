% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/2_1_nomogram.R
\name{as_nomogram}
\alias{as_nomogram}
\title{Construct nomogram ojects for high-dimensional Cox models}
\usage{
as_nomogram(
  object,
  x,
  time,
  event,
  pred.at = NULL,
  fun.at = NULL,
  funlabel = NULL
)
}
\arguments{
\item{object}{Model object fitted by `hdnom::fit_*()` functions.}

\item{x}{Matrix of training data used for fitting the model.}

\item{time}{Survival time. Must be of the same length with
the number of rows as \code{x}.}

\item{event}{Status indicator, normally 0 = alive, 1 = dead.
Must be of the same length with the number of rows as \code{x}.}

\item{pred.at}{Time point at which to plot nomogram prediction axis.}

\item{fun.at}{Function values to label on axis.}

\item{funlabel}{Label for \code{fun} axis.}
}
\description{
Construct nomograms ojects for high-dimensional Cox models
}
\note{
The nomogram visualizes the model under the automatically
selected "optimal" hyperparameters (e.g. lambda, alpha, gamma).
}
\examples{
data(smart)
x <- as.matrix(smart[, -c(1, 2)])
time <- smart$TEVENT
event <- smart$EVENT
y <- survival::Surv(time, event)

fit <- fit_lasso(x, y, nfolds = 5, rule = "lambda.1se", seed = 1001)

nom <- as_nomogram(
  fit, x, time, event, pred.at = 365 * 2,
  funlabel = "2-Year Overall Survival Probability"
)

print(nom)
plot(nom)
}
