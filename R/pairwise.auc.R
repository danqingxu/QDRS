#' @title Pairwise Area Under the ROC Curve
#' @description \code{pairwise.auc} calculates AUROC for pairwise levels of group contrast .
#' @param x A response vector.
#' @param g A grouping vector or factor.
#' @param format Format of result. The default is "l", long format. The other option is square format.
#' @import PRROC
#' @export
#' @return A data frame with group 1, group 2, and the AUROC with comparison between them.
#' @examples
#' \dontrun{
#' data(EHR)
#' x = c(rnorm(2500),runif(2500))
#' res = pairwise.auc(x, EHR$sample.group)
#' }
pairwise.auc <- function(x,g,format="l"){
  group = levels(g)
  if(format=="l"){
    AUC = NULL
    group1 = NULL
    group2 = NULL
    for (i in 1:(length(group)-1)){
      for (j in (i+1):length(group)){
        group1 = c(group1,group[i])
        group2 = c(group2,group[j])
        roc.sub = c(group[i],group[j])
        dat.sub = data.frame(x=x[g %in% roc.sub],g=factor(g[g %in% roc.sub],levels=roc.sub))
        roc.obj = roc.curve(scores.class0=dat.sub$x, weights.class0 = (dat.sub$g==roc.sub[2]))
        AUC = c(AUC,roc.obj$auc)
      }
    }
    pairwise = data.frame(group1=group1,group2=group2,auc=AUC)
  } else {
    pairwise = data.frame(matrix(rep(NA,length(group)*length(group)),nc=length(group)))
    colnames(pairwise) = rownames(pairwise) = group
    for (i in 1:(length(group)-1)){
      for (j in (i+1):length(group)){
        roc.sub = c(group[i],group[j])
        dat.sub = data.frame(x=x[g %in% roc.sub],g=factor(g[g %in% roc.sub],levels=roc.sub))
        roc.obj = roc.curve(scores.class0=dat.sub$x, weights.class0 = (dat.sub$g==roc.sub[2]))
        pairwise[j,i] = roc.obj$auc
      }
    }
  }
  return(pairwise)
}
