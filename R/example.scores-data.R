#' @title An example resulting score set
#' @description It contains \code{score.mat} a matrix of five types of QDRSs for 5,000 subjects, \code{score.names} a vector of score names, \code{group} a vector of group with two levels "Case" and "Control", and \code{group.levels} the unique grouping values.
#' @docType data
#' @usage data(example.scores)
#' @examples
#' data(example.scores)
#' \dontrun{
#' # check the number of cases and controls
#' table(example.scores$.group)
#' }
#' @keywords datasets
"example.scores"
