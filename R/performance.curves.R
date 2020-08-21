#' @title Multiple ROC curves
#' @description Generate a pdf file with multiple plots for several group level contrast. Each plot consists of ROC curves for different quantitative disease scores, and the corrsponding AUROCs for performance comparison.
#' @param disease A string of disease name.
#' @param output.date An output date. The default is system date (today).
#' @param score.mat A resulting score matrix with multiple types of scores.
#' @param score.names A vector that indicates the scores of interest for comparison.
#' @param group A grouping vector.
#' @param group.levels A vector of grouping levels that reflects the severity from low to high. The typical first element represents control or healthy status.
#' @param legend.pos A vector of legend position. The default is "bottomright" for all plots.
#' @param pairs.sub An index vector for subset of group level pairs, e.g. the pairs contrasting control vs. other case stages.
#' @export
#' @import PRROC
#' @return It produces a pdf file of plots.
#' @examples
#' \dontrun{
#' multiple.roc.curve(disease = "Disease",
#'   score.mat = example.scores$score.mat,
#'   score.names = example.scores$score.names,
#'   group = example.scores$group,
#'   group.levels = example.scores$group.levels,
#'   legend.pos = "bottomright", pairs.sub = 1:5)
#' }

multiple.roc.curve <- function(disease, output.date = NULL, score.mat, score.names, group, group.levels, legend.pos = "bottomright", pairs.sub = NULL){
  if (length(score.names) > 10){
    message("Showing more than 10 scores in one plot is not recommended. Only first 10 scores will be used.")
    score.names = score.names[1:10]
  }
  group.counts = table(group)
  check.counts = as.vector(table(group)<10)
  if (any(check.counts)){
    stop(paste0(names(group.counts)[which(check.counts)], " group(s) has/have less than 10 subjects."))
  }
  if(is.null(output.date)) {output.date = as.character(Sys.Date())}

  if (length(group.levels) == 2 ){
    group.pairs = data.frame(group1 = group.levels[1],
                             group2 = group.levels[2])
  } else {
    group1 = rep(group.levels[-length(group.levels)], times = (length(group.levels)-1):1)
    group2 = NULL
    for (d in 1:(length(group.levels)-1)){
      group2 = c(group2, group.levels[-(1:d)])
    }
    group.pairs = data.frame(group1 = group1, group2 = group2)
  }
  if (is.null(pairs.sub)) {pairs.sub = 1:nrow(group.pairs)}
  plot.num = length(pairs.sub)

  ynum = ceiling(plot.num/2)
  xnum = ifelse(plot.num == 1, 1, 2)
  pdf.width = ifelse(plot.num == 1,8.5/2,8.5)
  switch (ynum,
    "1" = {pdf.height <- 11/3*1.2},
    "2" = {pdf.height <- 11/3*2},
    "3" = {pdf.height <- 11/3*3}
  )
  group = as.character(group)
  score.mat = as.data.frame(score.mat)
  score.mat = score.mat[!is.na(group),]
  group = group[!is.na(group)]
  check.names = !score.names%in%colnames(score.mat)
  if (any(check.names)){
    stop(score.names[which(check.names)]," do(es) not exist in the provided score matrix.")
  }
  if (length(legend.pos) == 1) {
    legend.pos = rep(legend.pos, plot.num)
  } else if (length(legend.pos) != plot.num) {
    message("The legend position vector does not have the same length as the number of plots. Only the first element is used.")
    legend.pos = rep(legend.pos[1], plot.num)
  }

  if (output.date == "none") {
    pdf(file = paste0(disease,"_roc_curves.pdf",))
  } else {
    pdf(file = paste0(disease,"_roc_curves_",output.date,".pdf"), height = pdf.height, width = pdf.width)
  }
  par(mfrow = c(ynum,xnum))
  roc.plots <- list()
  for (i in 1:length(pairs.sub)){
    g1 = group.pairs[pairs.sub[i],"group1"]
    g2 = group.pairs[pairs.sub[i],"group2"]
    score.roc = list()
    auc = rep(NA, length(score.names))
    for (j in 1:length(score.names)){
      score.roc[[j]] <- tmp <- roc.curve(scores.class0 = score.mat[group == g2, score.names[j]], scores.class1 = score.mat[group == g1, score.names[j]], curve = T)
      auc[j] = tmp$auc
      if (j == 1){
        plot(score.roc[[j]], max.plot = TRUE, min.plot = TRUE, rand.plot = TRUE,fill.area = T, color = curve.palettes[j], auc.main = FALSE, main = paste0(disease, " ", g2," vs. ",g1, " ROC curves"), lwd = 1.2, lty = j)
      } else {
        plot(score.roc[[j]], add = TRUE, color = curve.palettes[j], lwd = 1.2, lty = j);
      }
    }
    legend(x = legend.pos[i],
           legend = paste0(score.names, " AUC: ", round(auc,3)),
           col = curve.palettes[1:length(score.names)],
           lty = 1:length(score.names),
           lwd = 1.2, cex = 0.75)
  }
  dev.off()
}

#' @title Multiple PR curves
#' @description Generate a pdf file with multiple plots for several group level contrast. Each plot consists of PR curves for different quantitative disease scores, and the corrsponding AUPRs for performance comparison.
#' @param disease A string of disease name.
#' @param output.date An output date. The default is system date (today).
#' @param score.mat A resulting score matrix with multiple types of scores.
#' @param score.names A vector that indicates the scores of interest for comparison.
#' @param group A grouping vector.
#' @param group.levels A vector of grouping levels that reflects the severity from low to high. The typical first element represents control or healthy status.
#' @param legend.pos A vector of legend position. The default is "bottomright" for all plots.
#' @param pairs.sub An index vector for subset of group level pairs, e.g. the pairs contrasting control vs. other case stages.
#' @export
#' @import PRROC
#' @return It produces a pdf file of plots.
#' @examples
#' \dontrun{
#' multiple.pr.curve(disease = "Disease",
#'   score.mat = example.scores$score.mat,
#'   score.names = example.scores$score.names,
#'   group = example.scores$group,
#'   group.levels = example.scores$group.levels,
#'   legend.pos = "topright", pairs.sub = 1:5)
#' }

multiple.pr.curve <- function(disease, output.date = NULL, score.mat, score.names, group, group.levels, legend.pos = "topright", pairs.sub = NULL){
  if (length(score.names) > 10){
    message("Showing more than 10 scores in one plot is not recommended. Only first 10 scores will be used.")
    score.names = score.names[1:10]
  }
  group.counts = table(group)
  check.counts = as.vector(table(group)<10)
  if (any(check.counts)){
    stop(paste0(names(group.counts)[which(check.counts)], " group(s) has/have less than 10 subjects."))
  }
  if(is.null(output.date)) {output.date = as.character(Sys.Date())}

  if (length(group.levels) == 2 ){
    group.pairs = data.frame(group1 = group.levels[1],
                             group2 = group.levels[2])
  } else {
    group1 = rep(group.levels[-length(group.levels)], times = (length(group.levels)-1):1)
    group2 = NULL
    for (d in 1:(length(group.levels)-1)){
      group2 = c(group2, group.levels[-(1:d)])
    }
    group.pairs = data.frame(group1 = group1, group2 = group2)
  }
  if (is.null(pairs.sub)) {pairs.sub = 1:nrow(group.pairs)}
  plot.num = length(pairs.sub)

  ynum = ceiling(plot.num/2)
  xnum = ifelse(plot.num == 1, 1, 2)
  pdf.width = ifelse(plot.num == 1,8.5/2,8.5)
  switch (ynum,
          "1" = {pdf.height <- 11/3*1.2},
          "2" = {pdf.height <- 11/3*2},
          "3" = {pdf.height <- 11/3*3}
  )
  group = as.character(group)
  score.mat = as.data.frame(score.mat)
  score.mat = score.mat[!is.na(group),]
  group = group[!is.na(group)]
  check.names = !score.names%in%colnames(score.mat)
  if (any(check.names)){
    stop(score.names[which(check.names)]," do(es) not exist in the provided score matrix.")
  }
  if (length(legend.pos) == 1) {
    legend.pos = rep(legend.pos, plot.num)
  } else if (length(legend.pos) != plot.num) {
    message("The legend position vector does not have the same length as the number of plots. Only the first element is used.")
    legend.pos = rep(legend.pos[1], plot.num)
  }

  if (output.date == "none") {
    pdf(file = paste0(disease,"_pr_curves.pdf",))
  } else {
    pdf(file = paste0(disease,"_pr_curves_",output.date,".pdf"), height = pdf.height, width = pdf.width)
  }
  par(mfrow = c(ynum,xnum))
  PR.plots <- list()
  for (i in 1:length(pairs.sub)){
    g1 = group.pairs[pairs.sub[i],"group1"]
    g2 = group.pairs[pairs.sub[i],"group2"]
    score.pr = list()
    auc = rep(NA, length(score.names))
    for (j in 1:length(score.names)){
      score.pr[[j]] <- tmp <- pr.curve(scores.class0 = score.mat[group == g2, score.names[j]], scores.class1 = score.mat[group == g1, score.names[j]], curve = T)
      auc[j] = tmp$auc.integral
      if (j == 1){
        plot(score.pr[[j]], max.plot = TRUE, min.plot = TRUE, rand.plot = TRUE,fill.area = T, color = curve.palettes[j], auc.main = FALSE, main = paste0(disease, " ", g2," vs. ",g1, " PR curves"), lwd = 1.2, lty = j)
      } else {
        plot(score.pr[[j]], add = TRUE, color = curve.palettes[j], lwd = 1.2, lty = j);
      }
    }
    legend(x = legend.pos[i],
           legend = paste0(score.names, " AUPR: ", round(auc,3)),
           col = curve.palettes[1:length(score.names)],
           lty = 1:length(score.names),
           lwd = 1.2, cex = 0.75)
  }
  dev.off()
}
