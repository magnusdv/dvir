% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dviJoint.R
\name{swapRoles}
\alias{swapRoles}
\title{Swap centricity of an assignment table}
\usage{
swapRoles(df, from = NULL, to = NULL)
}
\arguments{
\item{df}{A data frame. Each row is an assignment, with \code{*} representing
nonpairing.}

\item{from}{A character vector; either victims or missing persons. By
default, the column names of \code{df}. The only time this argument is needed,
if when \code{df} has other columns in addition, as in output tables of
\code{dviJoint()}.}

\item{to}{The column names of the transformed data frame. If missing, the
unique elements of \code{df} are used. An error is raised if \code{to} does not
contain all elements of \code{df} (except \code{*}).}
}
\value{
A data frame with \code{nrow(df)} rows and \code{length(to)} columns.
}
\description{
This function switches the roles of victims and missing persons in a table of
assignments, from PM-centric (victims as column names) to AM-centric (missing
persons as column names), and \emph{vice versa}. In both version, each row
describes the same assignment vector.
}
\examples{
df = example1 |> generatePairings() |> expand.grid.nodup()
df
swapRoles(df)

# Swap is idempotent
stopifnot(identical(swapRoles(swapRoles(df)), df))

}
