#' @title An example data set
#' @description It contains \code{sample.set} a set of 95 features for 5,000 subjects, \code{sample.group} a vector of group with two levels "Case" and "Control", and \code{training} a logical vector that indicates the usage for training.
#' @docType data
#' @usage data(EHR)
#' @examples
#' data(EHR)
#' \dontrun{
#' # check the number of cases and controls
#' table(EHR$sample.group)
#' # check the prevalence of features
#' colMeans(EHR$sample.set)
#' }
#' @keywords datasets
"EHR"
