#' @title Predict Quantitative Disease Risk Score
#' @description Predicts Quantitative Disease Risk Scores using weights derived from training set.
#' @param Y A new data set for prediction.
#' @param weights A matrix of weights.
#' @export
#' @return It returns a matrix of resulting quantitative disease risk scores.
#' @examples
#' \dontrun{
#' data(EHR)
#' res1 <- PC(X = EHR$sample.set,
#'   group = EHR$sample.group,
#'   training = EHR$training, pc.num = 1:2)
#' res2 <- predictQDRS(Y = EHR$sample.set[seq(1,50),], weights = res1$weights)
#' }


predictQDRS <- function(Y, weights){
  if (is.null(dim(weights))) { weights = matrix(weights, nc = 1)}
  if (ncol(Y) != nrow(weights)) { stop("The dimensions of Y and weights do not match.") }
  result = as.matrix(Y)%*%as.matrix(weights)
  return(result)
}
