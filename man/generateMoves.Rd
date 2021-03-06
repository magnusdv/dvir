% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generateMoves.R
\name{generateMoves}
\alias{generateMoves}
\title{Generate list of possible moves}
\usage{
generateMoves(from, to, ids.to)
}
\arguments{
\item{from}{A list of singletons.}

\item{to}{A list of pedigrees.}

\item{ids.to}{Character vector with names of missing persons.}
}
\value{
A list of moves
}
\description{
A a sex consistent list of moves for `global` is generated
}
\details{
Identity moves are included.
}
\examples{
\dontrun{
from = list(singleton("V1", 1), singleton("V2", 2))
ids.to = paste("MP", 1:4, sep = "")
to = list(nuclearPed(children = ids.to[1:3]),
          nuclearPed(children = ids.to[4], sex = 2))
generateMoves(from, to,  ids.to)
}
}
\seealso{
`global`.
}
