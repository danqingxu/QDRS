#' Compute Phenotype Risk Scores (PheRSs).
#' @description \code{PheRS} computes Phenotype Risk Scores based on a approach that was proposed in the context of rare Mendelian phenotypes.
#' @param X The original data set that include training and test sets. It should be a matrix of binary numbers.
#' @param feature.prevalence A vector of feature prevalence.
#' @param group A vector that indicate cases ("Case") and controls ("Control"). \code{NA} is allowed. The default is NULL, meaning all observations will be used in prevalence computation if the feature prevalence vector is not provided by the user. Otherwise, only the controls will be used.
#' @export
#' @return weights The PheRS weights for input features.
#' @return scores The resulting PheRSs for the whole set.
#' @references Lisa Bastarache, Jacob J Hughey, Scott Hebbring, Joy Marlo, Wanke Zhao, Wanting T Ho, Sara L Van Driest, Tracy L McGregor, Jonathan D Mosley, Quinn S Wells, et al. Phenotype risk scores identify patients with unrecognized mendelian disease patterns. Science, 359(6381):1233– 1239, 2018.
#' @examples
#' \dontrun{
#' data(EHR)
#' res1 <- PheRS(X = EHR$sample.set, group = EHR$sample.group)
#' res2 <- PheRS(X = EHR$sample.set, group = NULL)
#' }
PheRS <- function(X,feature.prevalence=NULL,group=NULL){
  if(is.null(feature.prevalence)){
    if (is.null(group)){
      feature.prevalence=colMeans(X==1)
    } else {
      controls = X[group=="Control",]
      feature.prevalence=colMeans(controls==1)
    }
  }
  wt = log(1/feature.prevalence)
  PheRS = (X==1)%*%log(1/feature.prevalence)
  return(list(weights=wt,score=PheRS))
}
