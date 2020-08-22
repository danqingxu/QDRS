#' @title Distribution Plots of QDRS
#' @description Generate a pdf file with a density plot, a boxplot and a prevalence plot for the resulting scores. These three plots show the distribution of scores by disease status for visualization of risk stratification.
#' @param disease A string of disease name.
#' @param output.date An output date. The default is system date (today).
#' @param score.mat A resulting score matrix with multiple types of scores.
#' @param score.name The score of interest.
#' @param group A grouping vector.
#' @param cutoff A vector of cutoff points to set the bins of prevalence plot. The default vector has 60 bins.
#' @param unknown.show A logical value indicates whethe unknowns should be used for plotting.
#' @export
#' @import PRROC
#' @return It produces a pdf file of plots and returns the list of three distribution plots.
#' @examples
#' \dontrun{
#' group1 = example.scores$group
#' group1[group1 != "Control"] = "Case"
#' set.seed(830)
#' na.ind = sample(1:length(group), size = 1000)
#' group1[na.ind] = NA
#' res = dist.plot(disease = "Disease",
#'   score.mat = example.scores$score.mat,
#'   score.name = "LPC",
#'   group = group1,
#'   unknown.show = T)
#' }

dist.plot <- function(disease, output.date = NULL, score.mat, score.name, group, cutoff = c(seq(0.05,0.45,by=0.05),seq(0.5,1,by=0.01)), unknown.show = FALSE){
  if (length(score.name) > 1){
    message("More than one score are specified. Only the first score is plotted.")
    score.name = score.name[1]
  }
  if (length(unique(group[!is.na(group)]))!=2){
    stop("group can only take on values 'Case', 'Control', or NA.")
  }
  group.counts = table(group)
  check.counts = as.vector(table(group)<30)
  if (any(check.counts)){
    stop(paste0(names(group.counts)[which(check.counts)], " group(s) has/have less than 30 subjects."))
  }
  if(is.null(output.date)) {output.date = as.character(Sys.Date())}

  unknown = FALSE
  if (unknown.show) {
    if (sum(is.na(group))==0) {
      message("No unknown/missing observations.")
    } else {
      group[is.na(group)] = "Unknown"
      unknown = TRUE
    }
  }

  prev_plots_set = data.frame(score.mat, Group = group, label = (group == "Case"))
  #---------------- a Density plots --------------#
  density_plot <- ggplot(prev_plots_set, aes_string(x = score.name)) +
    geom_density(color="darkblue", fill="lightblue") +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5),
          panel.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.ticks.x=element_blank(),
          text = element_text(size = 10)) +
    labs(x = paste0(score.name, " for ", disease), y="Density",title = "")
  #---------------- b Boxplots --------------#
  RS.Percentile = prev_plots_set[,c("Group","label")]
  RS.Percentile = cbind(RS.Percentile, ecdf(prev_plots_set[,score.name])(prev_plots_set[,score.name])*100)
  colnames(RS.Percentile)[3] = score.name
  if (unknown) {
    box_plot <- ggplot(subset(RS.Percentile,!is.na(Group)),aes_string(x = "Group",y = score.name, color = "Group", fill = "Group")) +
      geom_boxplot(outlier.shape = NA,outlier.size = 1) +
      theme_classic() +
      theme(legend.position = "none",
            plot.title = element_text(hjust = 0.5),
            panel.background = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.ticks.x=element_blank(),
            text = element_text(size=10)) +
      labs(x = disease, y = paste0(score.name, " percentile"), title = "") +
      scale_color_manual(values=c("black","black","black")) +
      scale_fill_manual(values=c("#2171b5","#6baed6","#c6dbef"))
  } else {
    box_plot <- ggplot(subset(RS.Percentile,!is.na(Group)),aes_string(x = "Group",y = score.name, color = "Group", fill = "Group")) +
      geom_boxplot(outlier.shape = NA,outlier.size = 1) +
      theme_classic() +
      theme(legend.position = "none",
            plot.title = element_text(hjust = 0.5),
            panel.background = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.ticks.x=element_blank(),
            text = element_text(size=10)) +
      labs(x = disease, y = paste0(score.name, " percentile"), title = "") +
      scale_color_manual(values=c("black","black")) +
      scale_fill_manual(values=c("#2171b5","#c6dbef"))
  }
  #---------------- c Prevalence Plot --------------#
  quantile.cutoff = data.frame(cutoff = c(0,cutoff), qt = quantile(prev_plots_set[,score.name], probs = c(0,cutoff)))
  preval = sum(prev_plots_set$label[prev_plots_set[,score.name]>=quantile.cutoff[1,"qt"]&prev_plots_set[,score.name]<=quantile.cutoff[2,"qt"]],na.rm=T)/sum(prev_plots_set[,score.name]>=quantile.cutoff[1,"qt"]&prev_plots_set[,score.name]<=quantile.cutoff[2,"qt"])
  for (i in 2:(nrow(quantile.cutoff)-1)){
    preval = c(preval,sum(prev_plots_set$label[prev_plots_set[,score.name]>=quantile.cutoff[i,"qt"]&prev_plots_set[,score.name]<=quantile.cutoff[i+1,"qt"]],na.rm=T)/sum(prev_plots_set[,score.name]>=quantile.cutoff[i,"qt"]&prev_plots_set[,score.name]<=quantile.cutoff[i+1,"qt"]))
  }
  percent = c(0,cutoff*100)
  perc.levels= paste(percent[-length(percent)],"-",percent[-1],sep="")
  phenotype.prevalence = data.frame(bins = perc.levels,
                                    percentile = cutoff*100,
                                    preval = preval)

 prevalence_plot <- ggplot(phenotype.prevalence, aes(x = percentile, y = preval*100, col = preval)) + geom_point() +
    theme_classic() +
    theme(legend.position="none",
          plot.title = element_text(hjust = 0.5),
          panel.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.ticks.x=element_blank(),
          text = element_text(size=10)) +
    labs(x = paste0("Percentile of ",score.name),y = paste0("Prevalence of ", disease," (%)"),title = "") +
    scale_color_gradient(low = "#4292c6", high = "#08306b")
 #---------------- Arrange --------------#
 figure <- ggarrange(density_plot, box_plot, prevalence_plot,
                     ncol = 3, nrow = 1)
 if (unknown) {grand.title = paste0(disease, ", with Unknowns")} else {
   grand.title = disease
 }
  figure <- annotate_figure(figure,
                 top = text_grob(grand.title, color = "#0033A0", face = "bold"))
 ggexport(figure, filename = paste0(score.name, "_", disease,"_distribution_plots_",output.date,".pdf"), height = 8/3, width = 8)
 return(list(density = density_plot, box = box_plot, prev = prevalence_plot))
}
