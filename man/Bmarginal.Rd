% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Bmarginal.R
\name{Bmarginal}
\alias{Bmarginal}
\title{Posterior pairing probabilities}
\usage{
Bmarginal(jointRes, missing, prior = NULL)
}
\arguments{
\item{jointRes}{Output from \code{\link[=jointDVI]{jointDVI()}}.}

\item{missing}{Character vector with names of missing persons.}

\item{prior}{A numeric vector of length equal the number of rows in
\code{jointRes}. Default is a flat prior.}
}
\value{
A matrix. Row \code{i} gives the posterior probability that victim \code{i} is
one of the missing persons or someone else, denoted '*'.
}
\description{
Compute posterior pairing and non-pairing probabilities, based on a prior and
the output from \code{\link[=jointDVI]{jointDVI()}}.
}
\details{
The prior assigns a probability to each assignment, each row of \code{jointRes}.
If the prior is not specified, a flat prior is used. The prior needs not sum
to 1 since the user may rather choose a flat prior on the \emph{a priori} possible
assignments.
}
\examples{
jointRes = jointDVI(example1)

Bmarginal(jointRes, example1$missing)

# Artificial example: all but optimal solution excluded by prior
Bmarginal(jointRes, example1$missing, prior = c(1, rep(0,26)))


}
\seealso{
\code{\link[=jointDVI]{jointDVI()}}
}
