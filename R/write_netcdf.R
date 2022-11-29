#' Function to write results from downscalr into NetCDF file format
#'
#'
#' @param res Result from downscale
#' @param rasterfile RasterLayer object with ns as values
#' @param filename Name and path of NetCDF file
#' @param raster.crs CRS of raster
#'
#' @return
#' * \file{write_netcdf} with name and path specified in filename function argument

#'
#' @export write_netcdf
#'
#' @importFrom raster extent area setValues getValues res as.matrix
#' @import ncdf4
#' @import sp
#'
#' @examples
#' require(dplyr)
#' require(tidyr)
#' require(tibble)
#' require(ncdf4)
#' require(raster)
#' 
#' ## A basic example to plot the observed LU changes
#'
#' areas <- data.frame(ns=getValues(argentina_raster),
#'                     area=getValues(area(argentina_raster))) %>%
#'          na.omit() %>%
#'          group_by(ns) %>%
#'          summarise(area=sum(area)) %>%
#'          mutate(ns=as.character(ns))
#'
#' to.plot <- argentina_luc %>%
#'            rename("times"="Ts") %>%
#'            left_join(areas) %>%
#'            mutate(value=value*area) %>%
#'            dplyr::select(-area)
#'
#' \dontrun{write_netcdf(to.plot ,argentina_raster, filename="argentina.nc")}
#'
write_netcdf <- function(res,
                         rasterfile,
                         filename = "nc_v1.nc",
                         raster.crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0") {
  # define NULL values because package error checking does not deal with dplyr correctly
  times = ns = lu.to = value = label = . = NULL
  
  p4s = sp::CRS(raster.crs)
  na_sum = function(x) {
    return(sum(x[!is.na(x)]))
  }
  geosims <- rasterfile
  sp::proj4string(geosims) <- p4s
  save_geovals <- data.frame(raster::getValues(geosims))
  colnames(save_geovals) <- "ns"
  
  geosims_area <- raster::area(geosims)
  ttemp = raster::getValues(geosims_area)
  ttemp[is.na(raster::getValues(geosims))] = NA
  geosims_area = raster::setValues(geosims_area, ttemp)
  
  
  
  # Define some straightforward dimensions
  lon <- seq(
    raster::extent(geosims)[1] +
      (raster::res(geosims)[2] / 2),
    raster::extent(geosims)[2] -
      (raster::res(geosims)[2] / 2),
    length.out = dim(geosims)[2]
  )
  lat <- seq(
    raster::extent(geosims)[4] - (raster::res(geosims)[1] / 2),
    raster::extent(geosims)[3] + (raster::res(geosims)[1] / 2),
    length.out = dim(geosims)[1]
  )
  time <- sort(unique(res$times))
  lc_class <- sort(unique(res$lu.from))
  
  
  dim_lon <- ncdf4::ncdim_def("lon", "degrees_east", lon)
  dim_lat <- ncdf4::ncdim_def("lat", "degrees_north", lat)
  dim_time <- ncdf4::ncdim_def("time", "years", time, unlim = TRUE)
  dim_lc_class <- ncdf4::ncdim_def("lc_class",
                                   paste(paste0(seq(
                                     1, length(lc_class), 1
                                   ), "=", lc_class), collapse = "/"), seq(1, length(lc_class), 1))
  dim_area <- ncdf4::ncdim_def("area", paste0("land_area"), 1)
  
  # create variables
  fillvalue <- NA
  LandCover_pixshare_long <-
    'share of pixel occupied by various land covers'
  LandCover_pixshare <-
    ncdf4::ncvar_def(
      name = "LC_area_share",
      units = "share of pixel area",
      longname = LandCover_pixshare_long,
      dim = list(dim_lon, dim_lat, dim_lc_class, dim_time),
      missval = fillvalue,
      prec = "double",
      compression = 9
    )
  
  PixelArea_long <- paste0('Total land area cover')
  PixelArea <-
    ncdf4::ncvar_def(
      name = paste0('Land_area'),
      units = "as input (km2)",
      longname = PixelArea_long,
      dim = list(dim_lon, dim_lat, dim_area),
      missval = fillvalue,
      prec = "double",
      compression = 9
    )
  
  
  
  ncid_out <-
    ncdf4::nc_create(filename, list(LandCover_pixshare, PixelArea))
  
  
  ### for the LU casef
  #add the declaration
  # reorder classes to match template
  res_file_LU <- res %>%
    group_by(times, ns, lu.to) %>%
    summarise(value = sum(value)) %>%
    rename(label = lu.to) %>%
    pivot_wider(
      id_cols = c(times, ns),
      names_from = "label",
      values_from = "value"
    ) %>%
    ungroup() %>%
    dplyr::mutate(across(!c(ns, times), ~ ifelse(is.na(.), 0, .))) %>%
    mutate(ns_area = rowSums(dplyr::select(., !c(ns, times)))) %>%
    mutate(across(!c(ns, times, ns_area), ~ . / ns_area))
  
  ns_area <- res_file_LU %>%
    subset(times == min(times)) %>%
    dplyr::select(c(ns, ns_area))
  
  res_file_LU <- res_file_LU %>%
    dplyr::select(-ns_area) %>%
    filter(times >= min(times)) %>%
    pivot_longer(
      cols = !c(ns, times) ,
      names_to = "label",
      values_to = "value"
    ) %>%
    group_by(times, ns, label) %>%
    arrange(label) %>%
    pivot_wider(
      id_cols = ns,
      names_from = c(times, label),
      values_from = value
    ) %>%
    left_join(ns_area) %>%
    mutate(ns = as.numeric(ns)) %>%
    ungroup()
  
  ns_area <- data.frame(ns_area)
  
  
  
  
  o = match(save_geovals$ns, res_file_LU$ns)
  to_plot2 = data.frame(res_file_LU[o, ])
  
  
  
  step = 2
  for (t in 1:length(time)) {
    for (c in 1:length(lc_class)) {
      temp_map <- raster::setValues(geosims, as.numeric(to_plot2[, step]))
      
      data <- raster::as.matrix(temp_map)
      ncdf4::ncvar_put(
        ncid_out,
        LandCover_pixshare,
        t(data),
        start = c(1, 1, c, t),
        count = c(ncol(data), nrow(data), 1, 1)
      )
      step = step + 1
    }
  }
  
  
  o2 = match(save_geovals$ns, res_file_LU$ns)
  to_plot_area = data.frame(ns_area[o2, ])
  temp_map <- raster::setValues(geosims, as.numeric(to_plot_area[, 2]))
  
  data <- raster::as.matrix(temp_map)
  
  ncdf4::ncvar_put(
    ncid_out,
    PixelArea,
    t(data),
    start = c(1, 1, 1),
    count = c(ncol(data), nrow(data), 1)
  )
  
  
  
  ncdf4::nc_close(ncid_out)
}
