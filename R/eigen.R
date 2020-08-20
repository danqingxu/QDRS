#' @title R Matrix with Rank One for Eigen Approach
#' @description Computes an estimated rank-one matrix R with unit-norm eigenvector v. Up to a sign ambiguity, the entries of v are proportional to the balanced accuracies (the average between the sensitivity and the specificity) of the input features.
#' @param Qmat The correlation matrix Q.
#' @import MASS
#' @return Rmat An estimate of the rank-one R matrix
rankOne.R <- function(Qmat) {
  L <- dim(Qmat)[1]
  N <- L*(L-1)/2

  # Set up system of equations
  A <- matrix( rep(0,L*N) ,ncol=L)
  q.vec <- rep(0,N)
  ind <- 1
  for (j in 1:(L-1)) {
    for (k in (j+1):L) {
      q.vec[ind] <- log(abs(Qmat[j,k]))
      A[ind,j] <- 1
      A[ind,k] <- 1
      ind <- ind + 1
    }
  }

  # Least squares solution
  A.inv <- ginv(A)
  t.vec <-  A.inv %*% q.vec
  # lm1 = lm(q.vec~A+0)
  # t.vec = as.vector(lm1$coefficients)

  # Calculate Rmat
  Rmat <- exp(t.vec) %*% t(exp(t.vec))
  return(Rmat)
}


#' @title Eigen Scores
#' @description Compute Eigen weights and scores. This method is unsupervised, and assigns weights that are proportional to the balanced accuracies (the average between the sensitivity and the specificity) of the input feature.
#' @param X The original data set that include training and test sets. It should be a matrix of numbers.
#' @param training A logical or index vector to indicate whether the subject belongs to the training set.
#' @return weights The Eigen weights for input features.
#' @return scores The resulting Eigen scores for the whole set.
#' @examples
#' \dontrun{
#' data(EHR)
#' res1 <- eigen.score(EHR$sample.set, EHR$training)
#' }
#' @references Iuliana Ionita-Laza, Kenneth McCallum, Bin Xu, and Joseph D Buxbaum. A spectral approach integrating functional genomic annotations for coding and noncoding variants. Nature genetics, 48(2):214, 2016.
eigen.score <- function(X, training){
  X = as.matrix(X)
  if(mode(X)=="character") { stop("The original data set must be numeric.") }

  Xmat <- scale(X)
  X.train <- Xmat[training,]
  cov1 <- cor(X.train,use="pairwise.complete.obs")
  Rmat <- rankOne.R(cov1)
  eigen.decomp <- eigen(Rmat)
  weights <- eigen.decomp$vectors[,1]
  wt = abs(weights/sum(weights))
  #wt = eigen.decomp$vectors[,1]
  score.N <- as.matrix(Xmat) %*% wt
  eigen.score <- list(weights = wt, scores = score.N)
  return(eigen.score)
}
