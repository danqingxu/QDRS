#' @title Pairwise Wilcoxon Rank Sum Tests
#' @description Calculates pairwise comparisons between group levels with corrections for multiple testing.
#' @param x A response vector.
#' @param g A grouping vector or factor.
#' @param alternative A character string specifying the alternative hypothesis, must be one of "two.sided" (default), "greater" or "less". You can specify just the initial letter.
#' @param p.adjust.method Method for adjusting p values (can be abbreviated): "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none". The default is "bonferroni".
#' @export
#' @return A data frame with group 1, group 2, and the Wilcoxon Rank Sum test p value of comparison between them.
#' @examples
#' \dontrun{
#' data(EHR)
#' x = c(rnorm(2500),runif(2500))
#' res = pairwise.wilcox(x, EHR$sample.group)
#' }

pairwise.wilcox <- function(x, g, alternative = "two.sided", p.adjust.method = "bonferroni"){
  wilcox.res = pairwise.wilcox.test(x, g, alternative = alternative, p.adjust.method = p.adjust.method)
  p = wilcox.res$p.value
  g1 = rep(as.character(colnames(p)),rep(nrow(p),ncol(p)))
  g2 = rep(as.character(rownames(p)),nrow(p))
  p = as.vector(p)
  res = data.frame(Group1=g1,Group2=g2,"p-value"=p)
  res = res[!is.na(p),]
  return(res)
}
