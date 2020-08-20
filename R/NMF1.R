#' @title Compute Non-negative Matrix Factorization with Rank 1 Scores (NMF1s).
#' @description \code{NMF1} can compute NMF1 scores with factor loadings constrained to be non-negative.
#' @param X The original data set. It should be a matrix of numbers.
#' @param p.seed A random seed for model fitting.
#' @export
#' @return lvs The resulting LVSs for the new set Y.
#' @examples
#' \dontrun{
#' data(EHR)
#' sample.set = EHR$sample.set
#' res1 <- NMF1(X = sample.set[seq(1,50),])
#' }
NMF1 <- function(X, p.seed = 123){
  rSums = rowSums(X)
  cSums = colSums(X)
  zero.rowind = which(rSums==0)
  zero.colind = which(cSums==0)
  message("Removing ", paste0(length(zero.rowind)," rows and ", length(zero.colind)," columns that have all zeros."))
  X.clean = X[-zero.rowind,-zero.colind]
  nmf1.out = nmf(X.clean, 1, seed = p.seed)
  scores = rep(0,nrow(X))
  scores[which(rSums!=0)] = as.numeric(nmf1.out@fit@W)
  weights = as.numeric(nmf1.out@fit@H)
  return(list(weights = weights, scores = scores))
}
