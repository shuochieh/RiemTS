#' CaliforniaRain
#' 
#' A real-data example of compositional time series. The data contain yearly
#' precipitation in California from 1895 to 2024, with each monthly value 
#' expressed as the percentage of total annual precipitation, i.e.,
#' precipitation in a given month divided by total precipitation for that year.
#' 
#' The (unnormalized) monthly precipitation is calculated as the mean 
#' precipitation observed at NOAA stations. Raw data are sourced from 
#' \url{https://www.ncei.noaa.gov/access/monitoring/climate-at-a-glance/divisional/time-series}
#' 
#' @examples 
#' data(CaliforniaRain)
#' str(CaliforniaRain)
"CaliforniaRain"
