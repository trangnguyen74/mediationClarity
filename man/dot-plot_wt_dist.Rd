% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Ax-shared_weights.R, R/A3-weights_odds.R,
%   R/C3a-estimate_psYpred.R
\name{.plot_wt_dist}
\alias{.plot_wt_dist}
\alias{.plot_wt_dist.odds}
\alias{.plot_wt_dist.psYpred}
\alias{dot-plot_wt_dist}
\title{(For maintainer) Plot weight distributions}
\usage{
.plot_wt_dist(w.dat, point.alpha = 0.1, jitter.width = 0.3)

.plot_wt_dist.odds(w.dat, point.alpha = 0.1, jitter.width = 0.3)

.plot_wt_dist.psYpred(w.dat, point.alpha = 0.1, jitter.width = 0.3)
}
\arguments{
\item{w.dat}{The weighted dataset.}

\item{point.alpha, jitter.width}{Graphical parameters for \code{geom_point()}.}
}
\value{
If constant sampling weights, plot the densities of distribution-morphing weights in the pseudo samples

If non-constant sampling weights, output two plots, one of the distribution-morphing weights, one of the final weights (product of sampling and distribution-morphing weights)
}
\description{
Internal function called by \code{weights_\*()} and \code{estimate_\*()} to make weight distribution plot(s).
}
\keyword{internal}
