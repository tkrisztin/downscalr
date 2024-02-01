#' Function to write results from downscalr or other sources into netCDF file format
#'
#'
#' @param data If not provided in \code{variables}, either a dataframe or matrix object with columns ns (raster cell numbers), times (optional), var1 (variable levels), and value; defaults to NULL
#' @param rasterfile Raster grid object with ns as cell values
#' @param filename Name of the output netCDF. Defaults to "Downscale_netCDF.nc"
#' @param variables A list of lists for each variable to be written or if \code{data} is provided a list of netCDF attribute settings as follows:
#' * \code{data} = a data frame object with columns ns (raster cell numbers), times (optional), var1 (variable levels), and value. Overridden if argument \code{data} is provided, defaults to NULL
#' * \code{name_long} = character string specifying a long descriptive name; defaults to name_long="This is long name"
#' * \code{name} = character string specifying the standard name, contains no white space and is case sensitive : defaults to name="this_is_standard_name"
#' * \code{units} = character string specifying the unit of measure of data; defaults to units="unit_of_measurement"
#' * \code{dimn} = integer specifying the number of dimensions; defaults to NULL
#' * \code{varn} = integer specifying the number of levels of variable; defaults to NULL
#' * \code{timen} = integer specifying the number of timesteps
#' * \code{datadim} = integer specifying the dimension of data; defaults to 1
#' * \code{expandValue} = numeric value or NA to fill in missing observations in \code{data}
#' * \code{vardescr} = character string providing longer description of variable levels in data; defaults to NULL
#' * \code{create.dimvar} = either TRUE or FALSE, defines if dimension variable of \code{data} should be created or not; defaults to TRUE
#' @param start.time Integer defining first time step of time dimension to write, if dimension \code{times} is provided in \code{data} or \code{variables}, ; defaults to min(data$times) per variable
#' @param end.time Integer defining last time step of time dimension to write, if dimension \code{times} is provided in \code{data} or \code{variables}, ; defaults to min(data$times) per variable
#' @param by.time Integer specifying interval of time steps to be written between \code{start.time} and \code{end.time}; defaults to 10
#' @param filename Character string of netCDF file name, has to include ".nc" extension; defaults to filename="Downscale_netCDF.nc"
#' @param filepath Character string to specify file path netCDF is written to, defaults to current working directory
#' @param verbose Either TRUE or FALSE, prints each layer written into the netCDF file (for troubleshooting purposes); defaults to FALSE
#' @return
#' * Outputs \file{netCDF} file with \code{filename} written into \code{filepath} as specified in function arguments
#'
#'
#' @export write_netcdf
#'
#'
#' @import ncdf4
#' @importFrom terra crs values ext res setValues
#' @import sp
#' @import dplyr
#'
#'
# @examples
#'
# start.year=2020
# label <- "label"
# units <- "km2"
# filename <- "./testMW.nc"
#
write_netcdf <- function(data=NULL, rasterfile=NULL, variables=list(name_long="This is long name", name="this_is_standard_name", units="unit_of_measurement", dimname="dim_name", dimn=3, varn=NULL, timen=NULL, datadim=1, expandValue=0, data=NULL, vardescr=NULL, create.dimvar=TRUE), start.time=NULL, end.time=NULL, by.time=10, filename = NULL, filepath=NULL, verbose=FALSE){

  times <- ns <- lu.to <- value <- var1 <-  NULL

  #input check
  #check whether variables input is list of list
  if(is.null(rasterfile)) stop(paste0("No raster provided!"))

  switch(all(sapply(variables,is.list)),
         switch(any(sapply(variables,function(x){is.null(x$data)})),stop(paste0("No data provided in one variable!"))),
         switch(is.null(data) | is.null(variables$data), stop(paste0("No data provided!")), NULL)
         )


   if(any(colnames(data) %in% c("times"))){
    if(is.null(start.time) & !is.null(data)) start.time <- min(data$times)
    if(is.null(end.time) & !is.null(data)) end.time <- max(data$times)
  }
  if(is.null(by.time)) by.time <- 10
  if(is.null(filename)) filename <- "Downscale_netCDF.nc"
  if(is.null(filepath)) filepath <- paste0(getwd())


  if(any(is.matrix(data), is.data.frame(data), !is.null(data), !is.null(variables$data))){
    filepath <- paste0(getwd())

    if (is.matrix(data)) data <- as.data.frame(data)
    if(is.list(variables)) data <- variables$data

    res <- data %>% dplyr::select(any_of(c("times","ns","value")),starts_with("var")) %>%     group_by(ns,across(starts_with("var")),times) %>% summarize(value=sum(value)) %>% ungroup()

    data_layer <- variables
    data_layer$data <- as.data.frame(res)
    data_layer$varn <- length(unique(res$var1))
    data_layer$data.dim <- length(grep("var",colnames(res)))

    if(any(colnames(res) %in% "times")) data_layer$timen <- length(unique(res$times))

    variables <- list(data_layer)

  }else{


  }

  on.exit(ncdf4::nc_close(ncid_out))


  p4s <- as.character(terra::crs(rasterfile, proj=TRUE))
  na_sum <- function(x) {return(sum(x[!is.na(x)]))}
  if(as.character(class(rasterfile))!="SpatRaster") geosims <- terra::rast(rasterfile) else geosims <- rasterfile
  #sp::proj4string(geosims) <- p4s
  if(p4s!="") geosims <- terra::project(geosims,p4s)
  save_geovals <- data.frame(terra::values(geosims))
  colnames(save_geovals) <- "ns"
  geosims_extent <- terra::ext(geosims)

  # Define some straightforward dimensions
  lon <- seq(geosims_extent[1]+(terra::res(geosims)[2]/2),geosims_extent[2]-(terra::res(geosims)[2]/2),length.out=dim(geosims)[2])
  lat <- seq(geosims_extent[4]-(terra::res(geosims)[1]/2),geosims_extent[3]+(terra::res(geosims)[1]/2),length.out=dim(geosims)[1])

  if(!is.null(start.time) & !is.null(end.time)){

    time <- sort(seq.int(start.time,end.time,by.time))
    dim_time <- ncdf4::ncdim_def("time", "years", calendar="standard",time,unlim=TRUE,longname=NULL)

                                        }


  dim_lon <- ncdf4::ncdim_def("longitude", "degrees_east",lon,longname=NULL)
  dim_lat <- ncdf4::ncdim_def("latitude", "degrees_north",lat,longname=NULL)




  fillvalue <- NA
  nc_var_list <- list()
  nc_plot_list <- list()
  label_name_list <- list()
  label_num_list <- list()

  for(ll in seq_len(length(variables))){


    layer_temp <- variables[[ll]]
    layer_temp$datadim <- length(grep("var", colnames(layer_temp$data)))
    if(length(layer_temp$dimname)<layer_temp$datadim & layer_temp$dimname=="dim_name") layer_temp$dimname <- paste0(layer_temp$dimname,seq_len(layer_temp$datadim))
    dim_temp_list <- list(dim_lon,dim_lat)
    label_names_temp_list <- list()

    for(kk in seq_len(layer_temp$datadim)){
      label_names_temp <- gsub("[.]"," ",unique(layer_temp$data[,paste0('var',kk)]))
      label_num_temp <- seq(1,length(label_names_temp),1)

      if(layer_temp$create.dimvar) dim_temp <- ncdf4::ncdim_def(paste0(layer_temp$dimname[kk],"_legend"),paste(paste(paste0(seq(1,length(label_names_temp),1)," = ",label_names_temp), collapse="; ")), label_num_temp, longname=NULL,create_dimvar=layer_temp$create.dimvar) else dim_temp <- ncdf4::ncdim_def(paste0(layer_temp$dimname[kk],"_legend"),"", as.integer(label_num_temp), longname=NULL,create_dimvar=FALSE)

      dim_temp_list[[kk+2]] <- dim_temp
      label_names_temp_list[[kk]] <- label_names_temp
    }
    # }else{
    #   kk <- 0
    #   }


    if(!any(is.na(layer_temp$timen),is.null(layer_temp$timen))){

      if(is.null(start.time)) start.time <- min(layer_temp$timen)
      if(is.null(end.time)) end.time <- max(layer_temp$timen)

      layer_temp_df <- layer_temp$data %>% filter(times>=start.time, times<=end.time) %>% ungroup() %>% group_by(times,across(starts_with("var"))) %>% arrange(times,across(starts_with("var")), .by_group=TRUE) %>%
        pivot_wider(id_cols = ns, names_from = c(times,grep('var',colnames(layer_temp$data),value=TRUE)), values_from = value)

      layer_temp_df[is.na(layer_temp_df)] <- layer_temp$expandValue

      layer_temp_df <- data.frame(layer_temp_df)


      dim_temp_list[[kk+3]] <- dim_time

      temp_nc_var <- ncdf4::ncvar_def(name = layer_temp$name,
                                      longname = layer_temp$name_long,
                                      dim = dim_temp_list, units = layer_temp$units,
                                      missval=fillvalue,prec = "double", compression = 9)


    }else{

      layer_temp_df <- layer_temp$data %>% ungroup() %>% group_by(across(all_of(grep('var',colnames(layer_temp$data),value=TRUE)))) %>% arrange(across(all_of(grep('var',colnames(layer_temp$data),value=TRUE))), .by_group=TRUE) %>%
        pivot_wider(id_cols = ns, names_from = c(grep('var',colnames(layer_temp$data),value=TRUE)), values_from = value)

      layer_temp_df[is.na(layer_temp_df)] <- layer_temp$expandValue

      layer_temp_df <- data.frame(layer_temp_df)


      temp_nc_var <- ncdf4::ncvar_def(name = layer_temp$name,units = layer_temp$units,
                                      longname = layer_temp$name_long,
                                      dim = dim_temp_list,
                                      missval=fillvalue,prec = "double",compression = 9)



    }
    # label and numbering
    o_temp = match(save_geovals$ns,as.integer(layer_temp_df[,c('ns')]))
    to_plot_temp = data.frame(layer_temp_df[o_temp,])
    colnames(to_plot_temp) <- gsub("[.]"," ",colnames(to_plot_temp))




    #define layer specific dim





    nc_var_list[[ll]] <- temp_nc_var
    nc_plot_list[[ll]] <- to_plot_temp
    label_name_list[[ll]] <- label_names_temp_list
    label_num_list[[ll]] <- label_num_temp
  }

  ncid_out <- ncdf4::nc_create(paste0(filepath,"/",filename),nc_var_list,force_v4=TRUE)

  ### add attributes

  for(pp in seq_len(length(nc_var_list))){
    nc_var_temp <- nc_var_list[[pp]]
    nc_plot_temp <- nc_plot_list[[pp]]
    #ncatt_put(ncid_out, nc_var_temp$name,"levels","gicht",prec="text")
    ncdf4::ncatt_put(ncid_out, nc_var_temp$name,"description",paste(paste(" ",variables[[pp]]$vardescr, collapse="; ")), prec="text")


    for (ii in 1:(ncol(nc_plot_temp)-1)) {
      t <- NULL ; #data_dim <- #c(rep(letters, )


      slct <- colnames(nc_plot_temp)[ii+1]
      if(!any(is.na(variables[[pp]]$timen), is.null(variables[[pp]]$timen))) t <- which(time==time[sapply(time,grepl,slct)]) else t <- NULL
      c <- unlist(lapply(label_name_list[[pp]],function(x){which(x==x[sapply(x,function(z) {grepl(z,slct)})])}))
      temp_map <- terra::setValues(geosims,nc_plot_list[[pp]][,slct])

      temp_data <- as.matrix(temp_map, wide=TRUE)
      #data <- as.matrix(temp_map)
      if(verbose){
        print(slct)
        print(temp_map)
        print(start_temp)
      }

      start_temp <- c(1,1,c,t)
      count_temp <- c(ncol(temp_data),nrow(temp_data), rep(1,length(start_temp)-2))
      ncdf4::ncvar_put(ncid_out, nc_var_temp, t(temp_data), start=start_temp,
                       count=count_temp)

      #print(step)

      #print(step)
    }



  }


  #ncatt_put(ncid_out,0, "title", path)
  ncdf4::ncatt_put(ncid_out,0, "institution", "IIASA")
  history <- paste("Leopold J.C. Ringwald", date(), sep=", ")
  ncdf4::ncatt_put(ncid_out,0, "history", history)
  ncdf4::ncatt_put(ncid_out,0, "raster metadata",  paste("CRS: ",p4s,"   Extent:", paste(geosims_extent[1:4],collapse=", "), " (xmin,xmax,ymin,ymax)", "   Resolution:", paste(terra::res(rasterfile), collapse=";")))


  #nc_close(ncid_out)

  #label_return <- list(label_table=label_tbl, label_map=label_full_map)
  #write.table(label_return$label_table, file=paste0("./outputs/LabelReadme_",gdx.path,'netCDF_LEO.txt'), quote=FALSE, sep="                                                       ", row.names=FALSE)
  #write.csv(label_return$label_map, file=paste0("./outputs/LabelMap_",gdx.path,'netCDF_LEO.csv'), row.names=FALSE)
}


