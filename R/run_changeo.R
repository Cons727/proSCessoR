#' Run Change-O In R
#'
#' The function runs the clonotyping analysis after successful generation of productive-T_sc-pass.tsv files, but failing to calculate threshold. Specifically, it calculates the distance to nearest neighbor, finds the threshold and groups clones. Saves the final output to heavy and light clone-pass.tsv to current directory. Immcantation libraries must be installed since Docker is not used.
#' @param heavy_path Path to heavy_productive-T_sc-pass.tsv.
#' @param light_path Path to light_productive-T_sc-pass.tsv.
#' @param project_name Name of the project.
#' @keywords clones
#' @export
#' @import alakazam
#' @import data.table
#' @import dplyr
#' @import ggplot2
#' @import scoper
#' @import shazam
#' @import airr
#'

run_changeo <- function(heavy_path, light_path, project_name){

  db <- bind_rows(
    lapply(c(heavy_path, light_path), airr::read_rearrangement),
    .id="input_id"
  )

  dist_nearest <- "dist_nearest"

  # Calculate distance-to-nearest
  db <- distToNearest(db,
                      sequenceColumn="junction",
                      vCallColumn="v_call",
                      jCallColumn="j_call",
                      cellIdColumn="cell_id",
                      model="ham",
                      first=FALSE,
                      normalize="len")
  output <- findThreshold(db$dist_nearest, method="gmm", model = "gamma-norm", cutoff = "user", spc = 0.995)
  title = "GMM"
  if(is.na(output@threshold)){
    output <- findThreshold(db$dist_nearest, method="density")
    title = "Density"
  }
  threshold <- output@threshold

  pdf("Processing_QC/threshold_plot.pdf")
  plot(output, title=title)
  dev.off()

  db <- hierarchicalClones(db,
                           threshold = threshold,
                           method = "nt",
                           cell_id = "cell_id",
                           clone = "clone_id",
                           only_heavy = FALSE,
                           split_light = TRUE,
                           summarize_clones = FALSE)


  writeChangeoDb(db %>% filter(input_id == 1) %>% select(-input_id), paste0("Processing_QC/", project_name, "_heavy_clone-pass.tsv"))
  writeChangeoDb(db %>% filter(input_id == 2) %>% select(-input_id), paste0("Processing_QC/", project_name, "_light_clone-pass.tsv"))

}
