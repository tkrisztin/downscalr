#' Plot function to visualize downscale results
#'
#'
#' @param res Result from downscale
#' @param rasterfile RasterLayer object
#' @param year Specify dates of results to plot (has to be in res)
#' @param LU Specify land-use of results to plot (has to be in res)
#' @param color Pixel color on the map
#' @param label Label of plot legend
#'
#' @return A list containing
#' * \code{LUC.plot} A ggplot object
#' * \code{plot.df} A dataframe used for constructing LUC.plot
#'
#' @export LUC_plot
#'
#' @import ggplot2
#' @import dplyr
#' @import ggthemes
#' @import methods
#'
#' @examples
#' ## A basic example
LUC_plot <- function(res, rasterfile, year=NULL, LU=NULL, color = "Greens", label = "Area in ha per pixel"){

  plot_df <- terra::as.data.frame(rasterfile, xy = TRUE, na.rm = FALSE)

  to.plot <- res$out.res

  colnames(plot_df) <- c("x", "y", "ns")

  ns = lu.to = times = value = x = y= NULL
  if(is.null(year) & is.null(LU)){
    inputs <- to.plot %>% dplyr::group_by(ns, lu.to, times) %>% dplyr::summarise(value = sum(value),.groups = "keep")
  } else if(!(is.null(year) | is.null(LU))){
    inputs <- to.plot %>% dplyr::group_by(ns, lu.to, times) %>% dplyr::summarise(value = sum(value),.groups = "keep") %>% subset(lu.to==LU & times==year)
  } else if(is.null(year)){
    inputs <- to.plot %>% dplyr::group_by(ns, lu.to, times) %>% dplyr::summarise(value = sum(value),.groups = "keep") %>% subset(lu.to==LU)
  } else {
    inputs <- to.plot %>% dplyr::group_by(ns, lu.to, times) %>% dplyr::summarise(value = sum(value),.groups = "keep") %>% subset(times==year)
  }

  plot_df <- merge(plot_df, inputs, by="ns")

  plot_obj <- ggplot2::ggplot() +
    ggplot2::geom_tile(data=plot_df, aes(x=x, y=y, fill=value, group=lu.to), alpha=0.8) +
    ggplot2::scale_fill_distiller(palette = color , name=label, direction = 1)+

    ggplot2::coord_equal() +
    ggthemes::theme_map() +
    ggplot2::theme(legend.position="bottom") +
    ggplot2::theme(legend.key.width=unit(2, "cm"))+
    if(is.null(year) & is.null(LU)){
      ggplot2::facet_grid(times~lu.to)
    } else if(!(is.null(year) | is.null(LU))){
      ;
    } else if(is.null(year)){
      ggplot2::facet_wrap(~times)
    } else {
      ggplot2::facet_wrap(~lu.to)
    }

  return(list(LUC.plot=plot_obj, plot.df=plot_df))
}
