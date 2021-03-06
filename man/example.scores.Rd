% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/example.scores-data.R
\docType{data}
\name{example.scores}
\alias{example.scores}
\title{An example resulting score set}
\format{
An object of class \code{list} of length 4.
}
\usage{
data(example.scores)
}
\description{
It contains \code{score.mat} a matrix of five types of QDRSs for 5,000 subjects, \code{score.names} a vector of score names, \code{group} a vector of group with two levels "Case" and "Control", and \code{group.levels} the unique grouping values.
}
\examples{
data(example.scores)
\dontrun{
# check the number of cases and controls
table(example.scores$.group)
}
}
\keyword{datasets}
