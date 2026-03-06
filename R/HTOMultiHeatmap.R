#' HTO MULTI Heatmap
#'
#' This function generates a heatmap of HTO expression per batch for method MULTI. It returns a ggplot object.
#' @param object Seurat object after HTO demultiplexing.
#' @param classification Identity for the heatmap (rownames of heatmap); for example "HTO_classification" or "HTO_name". Default "HTO_classification".
#' @param assay Assay to use data from; for PTID-Visit-SortPop, use "PTID_HTO". Default "HTO".
#' @keywords demux
#' @export
#'


HTOMultiHeatmap <- function(object, classification = "HTO_classification", assay = "HTO"){
  #classification <- "HTO_classification"
  ident.global <- "HTO_classification.global"
  object <- object
  batch <- levels(object@meta.data[["orig.ident"]])
  DefaultAssay(object) <- assay
  Idents(object = object) <- object[[classification, drop = TRUE]]
  classification <- object[[classification]]
  singlets <- which(x = object[[ident.global]] == 'Singlet')
  singlet.ids <- grep("Negative", sort(x = unique(x = as.character(x = classification[singlets, ]))),
                      value = TRUE, fixed = TRUE, invert = TRUE)
  doublets <- which(object[[ident.global]] == 'Multiplet')
  doublet.ids <- sort(x = unique(x = as.character(x = classification[doublets, ])))
  heatmap.levels <- c(singlet.ids, doublet.ids, 'Negative')
  object <- ScaleData(object = object, assay = assay, verbose = FALSE)
  data <- FetchData(object = object, vars = singlet.ids)
  Idents(object = object) <- factor(x = classification[, 1], levels = heatmap.levels)

  plot <- SingleRasterMap(
    data = data,
    raster = TRUE,
    feature.order = rev(x = singlet.ids),
    cell.order = names(x = sort(x = Idents(object = object))),
    group.by = Idents(object = object)
    ) +
    guides(color = FALSE) +
    ggtitle(batch)

  default.colors <- c("violetred", "darkblue", "lightblue")
  cols <- ifelse(object[[ident.global]] == "Multiplet", default.colors[2], ifelse(object[[ident.global]] == "Negative", default.colors[1], default.colors[3]))
  names(cols) <- colnames(object)
  cols <- cols[names(sort(Idents(object)))]

  # scale the height of the bar
  pbuild <- ggplot_build(plot = plot)
  y.range <- diff(x = pbuild$layout$panel_params[[1]]$y.range)
  y.pos <- max(pbuild$layout$panel_params[[1]]$y.range) + y.range * 0.015
  y.max <- y.pos + 0.02 * y.range
  x.min <- min(pbuild$layout$panel_params[[1]]$x.range) + 0.1
  x.max <- max(pbuild$layout$panel_params[[1]]$x.range) - 0.1
  plot <- plot +
    annotation_raster(
      raster = t(cols),
      xmin = x.min,
      xmax = x.max,
      ymin = y.pos,
      ymax = y.max
    ) +
    coord_cartesian(ylim = c(0, y.max), clip = 'off')


  plot <- plot + theme(line = element_blank(),
                       text = element_text(size = 15), axis.text.y = element_text(size = 23))
  return(plot)
}
