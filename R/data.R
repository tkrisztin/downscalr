#' Argentina land-use change driver data
#'
#' @name argentina_df
#' @docType data
#' @keywords argentina drivers, land use change
#'
#' Data of land-use change drivers, as well as observed land-uses in 2000 in Argentina
#'
#' @format A list with four items: \code{xmat}, \code{restrictions}, \code{pop_data}, \code{lu_levels}
#' \code{xmat} is a  data frame with 96400 rows and 3 variables:
#' \describe{
#'   \item{ns}{Pixel ids for Argentina}
#'   \item{ks}{Variable names}
#'   \item{value}{Value of variable}
#' }
#' \code{restrictions} is a  data frame with 38560 rows and 4 variables:
#' \describe{
#'   \item{ns}{Pixel ids for Argentina}
#'   \item{lu.from}{Origin land-use class of converion}
#'   \item{lu.to}{Destination land-use class to which lu.from is converted}
#'   \item{value}{Binary, 1 denotes that the conversion lu.from to lu.to is not allowed}
#' }
#' \code{pop_data} is a  data frame with 84832 rows and 4 variables, containing population projections until 2050
#' \describe{
#'   \item{ns}{Pixel ids for Argentina}
#'   \item{times}{Year of population projection}
#'   \item{ks}{totPop - total population or ruralPop - rural population}
#'   \item{value}{Number of persons}
#' }
#' \code{lu_levels} is a  data frame with 23136 rows and 3 variables, containing land-use observations in 2000
#' \describe{
#'   \item{ns}{Pixel ids for Argentina}
#'   \item{lu.from}{Land-use class}
#'   \item{value}{Ha of land}
#' }
NULL

#' FABLE Calculator projections for Argentina
#'
#' @name argentina_FABLE
#' @docType data
#' @keywords land use change projections
#'
#' Land-use change projections from the FABLE Calculator until 2050
#'
#' @format A data frame with 150 rows and 4 variables:
#' \describe{
#'   \item{times}{Time period of land-use changes}
#'   \item{lu.from}{Origin land-use class of converion}
#'   \item{lu.to}{Destination land-use class to which lu.from is converted}
#'   \item{value}{Ha land being converted from lu.from to lu.to in time period}
#' }
#' @source \url{https://www.abstract-landscapes.com/fable-calculator}
NULL

#' Land-use change in Argentina from 2000 to 2020
#'
#' @name argentina_luc
#' @docType data
#' @keywords landuse change
#'
#' Data of land-use changes in Argentina between 2000 and 2010, as well as 2010 and 2020.
#'
#' @format A data frame with 277632 rows and 5 variables:
#' \describe{
#'   \item{ns}{Pixel ids for Argentina}
#'   \item{lu.from}{Origin land-use class of converion}
#'   \item{lu.to}{Destination land-use class to which lu.from is converted}
#'   \item{Ts}{Time-period; either 2000 (changes 2000 to 2010) or 2010 (changes 2010 to 2020)}
#'   \item{value}{Percent of lu.from in pixel being converted}
#' }
#' @source MapBiomas
NULL

#' List of Argentina key raster objects
#'
#' @name argentina_raster
#' @docType data
#' @keywords argentina raster
#'
#' List of CRS object, raster values, extent object to reconstruct Argentina map.
#'
#' @format A list with 3 objects:
#' \describe{
#'   \item{CRS}{Coordinate Reference System object}
#'   \item{values}{A matrix with 431 rows and 237 cols, !NA elements grid location of Argentina}
#'   \item{extent}{Extent object for Argentina raster map}
#' }
#' @source GLOBIOM map
NULL
