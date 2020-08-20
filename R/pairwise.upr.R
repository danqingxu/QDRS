#' @title Pairwise Area Under the PR Curve
#' @description Calculates AUPRC for pairwise group levels contrast.
#' @param x A response vector.
#' @param g A grouping vector or factor.
#' @param format Format of result. The default is "l", long format. The other option is square format.
#' @import PRROC
#' @export
#' @return A data frame with group 1, group 2, and the AUPRC with comparison between them.
#' @examples
#' \dontrun{
#' data(EHR)
#' x = c(rnorm(2500),runif(2500))
#' res = pairwise.upr(x, EHR$sample.group)
#' }
pairwise.upr <- function(x,g,format="l"){
  group = levels(g)
  if(format=="l"){
    UPR = NULL
    group1 = NULL
    group2 = NULL
    for (i in 1:(length(group)-1)){
      for (j in (i+1):length(group)){
        group1 = c(group1,group[i])
        group2 = c(group2,group[j])
        pr.sub = c(group[i],group[j])
        dat.sub = data.frame(x=x[g %in% pr.sub],g=factor(g[g %in% pr.sub],levels=pr.sub))
        pr.obj = pr.curve(scores.class0=dat.sub$x, weights.class0 = (dat.sub$g==pr.sub[2]))
        UPR = c(UPR,pr.obj$auc.integral)
      }
    }
    pairwise = data.frame(group1=group1,group2=group2,upr=UPR)
  } else {
    pairwise = data.frame(matrix(rep(NA,length(group)*length(group)),nc=length(group)))
    colnames(pairwise) = rownames(pairwise) = group
    for (i in 1:(length(group)-1)){
      for (j in (i+1):length(group)){
        pr.sub = c(group[i],group[j])
        dat.sub = data.frame(x=x[g %in% pr.sub],g=factor(g[g %in% pr.sub],levels=pr.sub))
        pr.obj = pr.curve(scores.class0=dat.sub$x, weights.class0 = (dat.sub$g==pr.sub[2]))
        pairwise[j,i] = pr.obj$auc.integral
      }
    }
  }
  return(pairwise)
}
