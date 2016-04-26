#' @title Extracting from Multiple Raster Layers
#'
#' @keywords raster extract layers
#'
#' @param raster_vector A single raster layer or a vector
#' @param gis_data A data frame containing longitude and latitude pairs
#' @param projection User specified projection #may not do anything becasue extract may reset it
#' @param filename A user specified name that the data will be saved as, instead of the default raster_data.csv
#' @param storeInDirectory This boolean parameter asks wheter or not the user wants to store
#' the data in the working directory.
#'
#' @author Alton Barbehenn
#'
#' @return This function returns a data frame/matrix at each location specified by gis_data parameter.
#' If storeInDirectory is TRUE, then it also saves the data to a file (defaults to raster_data.csv).
#'
#' @export
#'
#' @import raster
#'
#' @note I'm using a raster brick because, "while a RasterBrick has to refer to one multi-layer file or
#' is in itself a multi-layer object with data loaded in memory, a RasterStack may 'virtually'
#' connect several raster objects written to different files or in memory. Processing will be
#' more efficient for a RasterBrick than for a RasterStack, but RasterStack has the advantage
#' of facilitating pixel based calculations on separate raster layers."
#' From http://geoscripting-wur.github.io/IntroToRaster/
#'
#' @examples ExtractMultiple(c(bio1, bio9), gis_data)
#' species_vector<-c("Abies amabilis", "Acer nigrum")
#' gis_data = BIEN.gis.species(species_vector)[2:3]
#' gis_data = na.omit(gis_data)
#' gis_data = unique(gis_data)
#' raster_vec <- c(raster("bio_1.bil"), raster("bio_9.bil"), raster("bio_4.bil"))
#' ExtractMultiple(raster_vector = raster_vec, gis_data = gis_data, filename = "your_file_name.csv", storeInDirectory = TRUE)
ExtractMultiple <- function(raster_vector, gis_data, projection = NULL, filename = NULL, storeInDirectory = FALSE) {
  #check to make sure raster_vector is not empty
  if(length(raster_vector) == 0) {
    stop("The raster* object passed into this function contains no raster layers.")
  }

  #it's not really an error if there is no gis data, it's just silly

  #define range
  lon_min = min(gis_data$longitude)
  lon_max = max(gis_data$longitude)
  lat_min = min(gis_data$latitude)
  lat_max = max(gis_data$latitude)

  #crop all vectors to speed everything up
  raster_vector <- lapply(raster_vector, crop, extent(c(lon_min, lon_max, lat_min, lat_max)))

  #make raster brick
  raster_brick <- brick(raster_vector)

  #DO I NEED TO MAKE PROJECTION A MORE GENERAL PARAMETER SO THAT I CAN HANDLE DIFFERENT PROJECTIONS PER MAP?
  #IT MAY NOT EVEN MATER BECAUSE THE EXTRACT MAY OVERRIDE IT ANYWAY...

  #set points to the correct projection using a default if no parameter is used
  if(!is.null(projection)) {
    str(gis_data)
    coordinates(gis_data) <- ~longitude+latitude
    proj4string(gis_data) <- CRS(projection)
  }
  else {
    str(gis_data)
    coordinates(gis_data) <- ~longitude+latitude
    proj4string(gis_data) <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
  }

  #extract from brick
  raster_data = extract(raster_brick, gis_data)

  #if user wants to store in the working directory
  if (storeInDirectory) {
    #write data to a file in the working directory
    if (!is.null(filename)) {
      #using user's chosen name
      write.csv(raster_data, file = filename)
    }
    else {
      #using default name
      write.csv(raster_data, file=paste0("raster_data.csv"))
    }
  }


  return(raster_data)
}
