#' @title Compute individual PCs (Principal Components) weights and scores with selected signs.
#' @description \code{PC} can compute individual PC scores and only requires weak labels to help select the signs of individual PCs.
#' @param X The whole original data set that include training and test sets. It should be a matrix of numbers.
#' @param group A vector that indicate cases ("Case") and controls ("Control"). \code{NA} is allowed, and means that the observation is not used in individual PC sign determination.
#' @param training A logical or index vector to indicate whether the subject belongs to the training set.
#' @param pc.num The vector of desired individual PCs.
#' @export
#' @return weights The selected PC weights for input features.
#' @return scores The resulting PC scores for the whole set.
#' @examples
#' \dontrun{
#' data(EHR)
#' res1 <- PC(X = EHR$sample.set, group = EHR$sample.group, training = EHR$training, pc.num = 1:2)
#' }

PC <- function(X, group, training, pc.num){
  X = as.matrix(X)
  if(mode(X)=="character") { stop("The original data set must be numeric.") }

  whole.mat = scale(X)
  training.mat = whole.mat[training,]
  training.group = group[training]
  cov1 = cor(training.mat, use="pairwise.complete.obs")
  eigen.res = eigen(cov1)
  lambdas = eigen.res$values
  u = eigen.res$vectors
  Pc = as.matrix(training.mat)%*%u
  pcs = data.frame(Pc)
  colnames(pcs) <- paste0("pc",1:(dim(pcs)[2]))

  lpc.n = max(pc.num)
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
  u.sub = uM[,pc.num,drop=F]
  PC = as.matrix(whole.mat)%*%u.sub
  scores = data.frame(PC)
  colnames(scores) <- paste0("pc",1:(dim(scores)[2]))
  return(list(weights=u.sub, scores=scores))
}
