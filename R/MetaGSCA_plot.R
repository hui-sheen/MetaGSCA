#' Draw Forest plot based on metaprop() result
#'
#' Draw Forest plot based on meta result and meta-p vectors on bootstrap runs
#'
#' @param meta output object of metaprop function on original set of datasets
#' @param meta.p p-value vector obtained by metaprop on bootstrap series of datasets
#' @param distinct.label bottom text line on forestplot with nperm, nboot, and p summary 
MetaGSCA_plot <- function(meta,meta.p,distinct.label) {
	forest(
		meta,
		leftlabs = c("Datasets", "No. of sig", "NPerm"),
		rightlabs = c("P statistic", "95%-CI", "Weight"),
		xlim = c(0, max(meta$upper)+0.1),
		digits = 4,
		col.diamond = "red",
		col.diamond.lines = "black",
		print.tau2 = FALSE,
		print.I2 = FALSE,
		prediction = FALSE,
		col.predict = "black",
		smlab = ""
	)
	grid::grid.text(label = distinct.label,
						x = grid::unit(0.12, "npc"), y =grid::unit(0.05, "npc"),
						just = "left", gp = grid::gpar(fontsize = 12, fontface = "bold"))
	grid::grid.text(label = paste(" ",format(round(median(meta.p, na.rm = TRUE),4), nsmall=4), " [",
													format(round(quantile(meta.p, probs=0.025, na.rm = TRUE),4), nsmall=4), ", ",
													format(round(quantile(meta.p, probs=0.975, na.rm = TRUE),4), nsmall=4), "]", sep=""),
						x = grid::unit(0.70, "npc"), y = grid::unit(0.05, "npc"),
						just = "left", gp = grid::gpar(fontsize = 12, fontface = "bold"))
}
