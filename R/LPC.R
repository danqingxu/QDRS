#' @title Linear Combination of Principal Components
#' @description Computes LPC weights and scores based on an almost unsupervised method that combines multiple PCs and only requires weak labels to help select the signs of individual PCs.
#' @param X The original data set that include training and test sets. It should be a matrix of numbers.
#' @param group A vector that indicate cases ("Case") and controls ("Control"). \code{NA} is allowed, and means that the observation is not used in individual PC sign determination.
#' @param training A logical or index vector to indicate whether the subject belongs to the training set.
#' @param scale A logical value to indicate whether the input features need to be scaled.
#' @export
#' @import AssocTests
#' @return It returns a list of following components:
#' \item{lpc.n}{The number of significant PCs (eigenvalues) suggested by Tracy-Widom test, a vector of the sign of individual PC.}
#' \item{pc.sign}{A vector of the sign of individual PC.}
#' \item{weights}{The LPC weights for input features.}
#' \item{scores}{The resulting LPC scores for the whole set.}
#' @examples
#' \dontrun{
#' data(EHR)
#' res1 <- LPC(X = EHR$sample.set, group = EHR$sample.group, training = EHR$training)
#' }

LPC <- function(X, group, training, scale = TRUE){
  X = as.matrix(X)
  if(mode(X)=="character") { stop("The original data set must be numeric.") }

  if (scale) {
    Xmat = scale(X)
  } else {
    Xmat = X
  }
  training.mat = Xmat[training,]
  training.group = group[training]
  cov1 = cor(training.mat, use="pairwise.complete.obs")
  eigen.res1 = eigen(cov1)
  lambdas = eigen.res1$values
  u = eigen.res1$vectors
  Pc = as.matrix(training.mat)%*%u
  pcs = data.frame(Pc)
  colnames(pcs) <- paste0("pc",1:(dim(pcs)[2]))
  tw.res1 = tw(eigenvalues = lambdas,eigenL = length(lambdas),criticalpoint = 0.9793)
  lpc.n = tw.res1$SigntEigenL
  if(lpc.n==0) { stop("Zero significant eigenvalues.") }

  if(is.null(training.group)){}
  pcs$Group = training.group

  pc.sign = rep(NA,lpc.n)
  pcsM = pcs[,c(1:lpc.n,ncol(pcs))]
  uM = u[,1:lpc.n,drop=F]
  for(i in 1:lpc.n){
    pc.sign[i] = sign(mean(pcs[which(pcs$Group=="Case"),i])-mean(pcs[which(pcs$Group=="Control"),i]))
    pcsM[,i] = pc.sign[i]*pcsM[,i]
    uM[,i] = pc.sign[i]*uM[,i]
  }
  pc.num = 1:lpc.n
  pc.sub = pcsM[,pc.num]
  wt = lambdas[pc.num]
  #LPC.training = as.matrix(pc.sub)%*%matrix(wt,nc=1)
  u.sub = uM[,pc.num,drop=F]
  wts = u.sub%*%matrix(wt,nc=1)
  LPC = Xmat%*%as.matrix(wts)
  res = list(lpc.n = lpc.n, pc.sign = pc.sign, weights = wts, scores = LPC)
  return(res)
}
