% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mu.mnl.R
\name{sqr_diff.mnl}
\alias{sqr_diff.mnl}
\title{Squared differences optimization function with multinomial logit}
\usage{
sqr_diff.mnl(x, mu, areas, targets, restrictions = NULL, cutoff = 0)
}
\arguments{
\item{x}{Bias correction parameters}

\item{mu}{Means of each class}

\item{areas}{Vector of areas}

\item{targets}{Targets to match}

\item{restrictions}{Optional restrictions. Binary, if 1, the log-odds are set to zero}

\item{cutoff}{Optional, defaults to zero. If set, all log-odds below cutoff are set to zero}
}
\value{
Squared differences of optimized values and targets
}
\description{
Squared differences optimization function with multinomial logit
}
\keyword{internal}
