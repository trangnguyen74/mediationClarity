% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/A1-weights_med.R, R/A2-weights_ipw.R,
%   R/A3-weights_odds.R, R/C3a-estimate_psYpred.R
\name{.get_smd.med}
\alias{.get_smd.med}
\alias{.get_smd.ipw}
\alias{.get_smd.odds}
\alias{.get_smd.psYpred}
\alias{dot-get_smd}
\title{(For maintainer) Compute  (standardized) mean differences}
\usage{
.get_smd.med(w.dat, c.vars, m.vars, c.std, m.std)

.get_smd.ipw(w.dat, vars, standardize)

.get_smd.odds(w.dat, vars, standardize)

.get_smd.psYpred(w.dat, c.vars, c.std)
}
\arguments{
\item{w.dat}{Dataset for pseudo samples, already cleaned and dummy coded.}

\item{c.vars}{Names of covariates, checked and dummied.}

\item{m.vars}{Names of mediators, checked and dummied.}

\item{c.std}{Covariates to be standardized, already checked.}

\item{m.std}{Mediators to be standardized, already checked.}

\item{cm.vars}{Names of covariates and mediators, checked and dummied.}

\item{cm.vars.std}{Covariates and mediators to be standardized, already checked.}
}
\value{
A data frame containing the (standardized) mean differences.
}
\description{
Internal functions called by \code{.plot_balance.} functions to compute (standardized) mean differences to be plotted.
}
\keyword{internal}
