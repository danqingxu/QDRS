% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PheRS.R
\name{PheRS}
\alias{PheRS}
\title{Phenotype Risk Score}
\usage{
PheRS(X, feature.prevalence = NULL, group = NULL)
}
\arguments{
\item{X}{The original data set that include training and test sets. It should be a matrix of binary numbers.}

\item{feature.prevalence}{A vector of feature prevalence.}

\item{group}{A vector that indicate cases ("Case") and controls ("Control"). \code{NA} is allowed. The default is NULL, meaning all observations will be used in prevalence computation if the feature prevalence vector is not provided by the user. Otherwise, only the controls will be used.}
}
\value{
It returns a list of following components:
\item{weights}{The PheRS weights for input features.}
\item{scores}{The resulting PheRSs for the whole set.}
}
\description{
Computes Phenotype Risk Scores based on a approach that was proposed in the context of rare Mendelian phenotypes.
}
\examples{
\dontrun{
data(EHR)
res1 <- PheRS(X = EHR$sample.set, group = EHR$sample.group)
res2 <- PheRS(X = EHR$sample.set, group = NULL)
}
}
\references{
Lisa Bastarache, Jacob J Hughey, Scott Hebbring, Joy Marlo, Wanke Zhao, Wanting T Ho, Sara L Van Driest, Tracy L McGregor, Jonathan D Mosley, Quinn S Wells, et al. Phenotype risk scores identify patients with unrecognized mendelian disease patterns. Science, 359(6381):1233– 1239, 2018.
}
