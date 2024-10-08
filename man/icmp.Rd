% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{icmp}
\alias{icmp}
\title{DVI dataset: A large reference pedigree}
\format{
A \code{dviData} object with the following content:
\itemize{
\item \code{pm}: A list of 5 singletons (victims).
\item \code{am}: A reference pedigree with 6 genotyped members and 12 missing
persons.
\item \code{missing}: A vector containing the names of the missing persons.
}
}
\usage{
icmp
}
\description{
DVI dataset based loosely on the ICMP 2017 workshop material
https://www.few.vu.nl/~ksn560/Block-III-PartI-KS-ISFG2017.pdf (page 18).
There are 3 female victims, 2 male victims and 6 missing persons of both
sexes. We have renamed the individuals and simulated data for 13 CODIS
markers (see Details).
}
\details{
The 13 markers are, in order: \code{CSF1PO}, \code{D3S1358}, \code{D5S818},\code{D7S820},
\code{D8S1179}, \code{D13S317}, \code{D16S539}, \code{D18S51}, \code{D21S11}, \code{FGA}, \code{TH01}, \code{TPOX},
and \code{vWA}.

Source code for the simulation, and a file containing the allele frequencies,
can be found in the \code{data-raw} folder of the GitHub repository:
https://github.com/magnusdv/dvir.
}
\examples{
icmp

# plotDVI(icmp)

# Markers and allele frequencies
db = pedtools::getFreqDatabase(icmp$pm)
db

}
\keyword{datasets}
