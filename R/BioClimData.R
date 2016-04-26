#' A Bio Clim data for a vector of species
#'
#' This function is utilizes functions from the BIEN3 Project Written by Brian Maitner 2015.
#' It takes as input a vector containing a set of species and a zone number. This function
#' then queries the BIEN3 database for all reported locations of a memeber of the species.
#' Then it creates a csv file containing the bioclim values of all these locations.
#'
#' As an optional parameters, the user can choose to enter thier own file name for the csv file,
#' or let the code use the default name. Also if you wish to use a .bil file in a different
#' the user can enter this as a string. Bio Clim data must be in the working directory or in the
#' directoery specified by the directory parameter.
#'
#' @param species_vector A vector of species in the BIEN3 database. This should be a vector of strings, for example species_vector<-c("Abies amabilis", "Acer nigrum").
#' @param zone_number The zone number is an integer 1-18 (inclusive) representing that specific Bio Clim, for example, Bio Clim 4 would be zone_number = 4.
#' @param filename The location of the exported data, if you choose to use that parameter. Defaults to NULL.
#' @param directory Location of the Bio Clim data. Defaults to NULL.
#' @param storeInDirectory A boolean parameter representing whether or not the exported data should be in the working directory (FALSE) or with the Bio Clim data (TRUE). Defaults to FALSE.
#'
#' @return This function returns the value(s) of a given Bio Clim raster layer at every reported point for a species vector.
#' This function exports a csv containing the same information that is returned.
#'
#' @export
#'
#' @note The examples are labeled as "Not run" becuase the function accesses files and if the directory
#' is not correct, when the code is tested it wont be able to find the .bil file.
#'
#' @examples
#' \dontrun{BIEN.BioClim.Species(species_vector = s_vec, zone_number = 7)
#' species_data <- BIEN.BioClim.Species(species_vector = s_vec, zone_number = 7, filename = "test.csv", directory = "/filelocation/", storeInDirectory = TRUE)}
#'
BIEN.BioClim.Species <- function(species_vector, zone_number, filename = NULL, directory = NULL, storeInDirectory = FALSE) {
  #invoke required packages
  require(sp)
  require(raster)
  require(stats)
  require(utils)

  #set working directory to the optional parameter if it's entered
  old_directory = getwd()
  if (!is.null(directory)) {
    setwd(directory)
  }

  #check to see if the BioClim file exists
  if (!file.exists(paste0("bio_", zone_number, ".bil"))) {
    stop("The file .bil file corrosponding to the zone number entered
         either does not exist or cannot be read.")
  }

  #take the input family vector and optimize it, then query BIEN3 for the known locations
  species_vector = unique(species_vector)
  species_location_data <- BIEN.gis.species(species_vector)

  #take the longitude and latitude data and clean it
  gis_data = species_location_data[2:3]
  gis_data = na.omit(gis_data)
  gis_data = unique(gis_data)

  #define bound region of raster map to speed up data processing
  lon_min = min(gis_data$longitude)
  lon_max = max(gis_data$longitude)
  lat_min = min(gis_data$latitude)
  lat_max = max(gis_data$latitude)

  #get the requested BioClim data
  bio_clim_data = raster(paste0("bio_", zone_number, ".bil"))
  #check to see if the BioClim data has been correctly extracted
  #fromDisk(bio_clim_data)

  #crop BioClim map to make it easier to precess
  bio_clim_data = crop(bio_clim_data, extent(c(lon_min, lon_max, lat_min, lat_max)))

  #correct data: it is stored as a decimal ten times larger than actual to maintain
  #precision and to save memory
  bio_clim_data <- bio_clim_data/10

  #make the position data have the correct longitude and latitude format for
  #the projection used by the raster data
  #str(gis_data)
  coordinates(gis_data) <- ~longitude+latitude
  #"The data is in the latitude/longitude coordinate reference system (not projected) and the datum is WGS84"
  #proj4string(gis_data) <- CRS("+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs")

  #extract BioClim data at each reported locations
  bio9dat=extract(bio_clim_data, gis_data)

  #if user want to store in the entered directory
  if (storeInDirectory) {
    #output
    if (!is.null(filename)) {
      write.csv(bio9dat, file = filename)
    }
    else {
      write.csv(bio9dat, file=paste0("bio_", zone_number, "_data.csv"))
    }
  }
  #otherwise store in old directory
  else {
    setwd(old_directory)
    #output
    if (!is.null(filename)) {
      write.csv(bio9dat, file = filename)
    }
    else {
      write.csv(bio9dat, file=paste0("bio_", zone_number, "_data.csv"))
    }
  }

  #restore working directory to original directory regardless of save preference and return data
  setwd(old_directory)
  return(bio9dat)
}



#' A Bio Clim data for a vector of genuses
#'
#' This function is utilizes functions from the BIEN3 Project Written by Brian Maitner 2015.
#' It takes as input a vector containing a set of genuses and a zone number. This function
#' then queries the BIEN3 database for all reported locations of a memeber of the genus
#' Then it creates a csv file containing the bioclim values of all these locations.
#'
#' As an optional parameters, the user can choose to enter thier own file name for the csv file,
#' or let the code use the default name. Also if you wish to use a .bil file in a different
#' the user can enter this as a string. Bio Clim data must be in the working directory or in the
#' directoery specified by the directory parameter.
#'
#' @param genus_vector A vector of species in the BIEN3 database. This should be a vector of strings, for example genus_vector<-c("Abutilon","Abronia").
#' @param zone_number The zone number is an integer 1-18 (inclusive) representing that specific Bio Clim, for example, Bio Clim 4 would be zone_number = 4.
#' @param filename The location of the exported data, if you choose to use that parameter. Defaults to NULL.
#' @param directory Location of the Bio Clim data. Defaults to NULL.
#' @param storeInDirectory A boolean parameter representing whether or not the exported data should be in the working directory (FALSE) or with the Bio Clim data (TRUE). Defaults to FALSE.
#'
#' @return This function returns the value(s) of a given Bio Clim raster layer at every reported point for a genus vector.
#' This function exports a csv containing the same information that is returned.
#'
#' @export
#'
#' @note The examples are labeled as "Not run" becuase the function accesses files and if the directory
#' is not correct, when the code is tested it wont be able to find the .bil file.
#'
#' @examples
#' \dontrun{BIEN.BioClim.Genus(species_vector = g_vec, zone_number = 7)
#' genus_data <- BIEN.BioClim.Genus(genus_vector = g_vec, zone_number = 7, filename = "test.csv", directory = "/filelocation/", storeInDirectory = TRUE)}
#'
BIEN.BioClim.Genus <- function(genus_vector, zone_number, filename = NULL, directory = NULL, storeInDirectory = FALSE) {
  #invoke required packages
  require(sp)
  require(raster)
  require(stats)
  require(utils)

  #set working directory to the optional parameter if it's entered
  old_directory = getwd()
  if (!is.null(directory)) {
    setwd(directory)
  }

  #check to see if the BioClim file exists
  if (!file.exists(paste0("bio_", zone_number, ".bil"))) {
    stop("The file .bil file corrosponding to the zone number entered
         either does not exist or cannot be read.")
  }

  #take the input family vector and optimize it, then query BIEN3 for the known locations
  genus_vector = unique(genus_vector)
  genus_location_data <- BIEN.gis.genus(genus_vector)

  #take the longitude and latitude data and clean it
  gis_data = genus_location_data[3:4]
  gis_data = na.omit(gis_data)
  gis_data = unique(gis_data)

  #define bound region of raster map to speed up data processing
  lon_min = min(gis_data$longitude)
  lon_max = max(gis_data$longitude)
  lat_min = min(gis_data$latitude)
  lat_max = max(gis_data$latitude)

  #get the requested BioClim data
  bio_clim_data = raster(paste0("bio_", zone_number, ".bil"))
  #check to see if the BioClim data has been correctly extracted
  #fromDisk(bio_clim_data)

  #crop BioClim map to make it easier to precess
  bio_clim_data = crop(bio_clim_data, extent(c(lon_min, lon_max, lat_min, lat_max)))

  #correct data: it is stored as a decimal ten times larger than actual to maintain
  #precision and to save memory
  bio_clim_data <- bio_clim_data/10

  #make the position data have the correct longitude and latitude format for
  #the projection used by the raster data
  #str(gis_data)
  coordinates(gis_data) <- ~longitude+latitude
  #"The data is in the latitude/longitude coordinate reference system (not projected) and the datum is WGS84"
  #proj4string(gis_data) <- CRS("+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs")

  #extract BioClim data at each reported locations
  bio9dat=extract(bio_clim_data, gis_data)

  #if user want to store in the entered directory
  if (storeInDirectory) {
    #output
    if (!is.null(filename)) {
      write.csv(bio9dat, file = filename)
    }
    else {
      write.csv(bio9dat, file=paste0("bio_", zone_number, "_data.csv"))
    }
  }
  #otherwise store in old directory
  else {
    setwd(old_directory)
    #output
    if (!is.null(filename)) {
      write.csv(bio9dat, file = filename)
    }
    else {
      write.csv(bio9dat, file=paste0("bio_", zone_number, "_data.csv"))
    }
  }

  #restore working directory to original directory regardless of save preference and return data
  setwd(old_directory)
  return(bio9dat)
}





#' A Bio Clim data for a vector of families
#'
#' This function is utilizes functions from the BIEN3 Project Written by Brian Maitner 2015.
#' It takes as input a vector containing a set of genuses and a zone number. This function
#' then queries the BIEN3 database for all reported locations of a memeber of the family
#' Then it creates a csv file containing the bioclim values of all these locations.
#'
#' As an optional parameters, the user can choose to enter thier own file name for the csv file,
#' or let the code use the default name. Also if you wish to use a .bil file in a different
#' the user can enter this as a string. Bio Clim data must be in the working directory or in the
#' directoery specified by the directory parameter.
#'
#' @param family_vector A vector of species in the BIEN3 database. This should be a vector of strings, for example family_vector<-c("Theaceae","Ericaceae").
#' @param zone_number The zone number is an integer 1-18 (inclusive) representing that specific Bio Clim, for example, Bio Clim 4 would be zone_number = 4.
#' @param filename The location of the exported data, if you choose to use that parameter. Defaults to NULL.
#' @param directory Location of the Bio Clim data. Defaults to NULL.
#' @param storeInDirectory A boolean parameter representing whether or not the exported data should be in the working directory (FALSE) or with the Bio Clim data (TRUE). Defaults to FALSE.
#'
#' @return This function returns the value(s) of a given Bio Clim raster layer at every reported point for a family vector.
#' This function exports a csv containing the same information that is returned.
#'
#' @export
#'
#' @note The examples are labeled as "Not run" becuase the function accesses files and if the directory
#' is not correct, when the code is tested it wont be able to find the .bil file.
#'
#' @examples
#' \dontrun{
#' BIEN.BioClim.Family(family_vector = f_vec, zone_number = 7)
#' genus_data <- BIEN.BioClim.Family(family_vector = f_vec, zone_number = 7, filename = "test.csv", directory = "/filelocation/", storeInDirectory = TRUE)
#' }
BIEN.BioClim.Family <- function(family_vector, zone_number, filename = NULL, directory = NULL, storeInDirectory = FALSE) {
  #invoke required packages
  require(sp)
  require(raster)
  require(stats)
  require(utils)

  #set working directory to the optional parameter if it's entered
  old_directory = getwd()
  if (!is.null(directory)) {
    setwd(directory)
  }

  #check to see if the BioClim file exists
  if (!file.exists(paste0("bio_", zone_number, ".bil"))) {
    stop("The file .bil file corrosponding to the zone number entered
         either does not exist or cannot be read.")
  }

  #take the input family vector and optimize it, then query BIEN3 for the known locations
  family_vector = unique(family_vector)
  family_location_data <- BIEN.gis.family(family_vector)

  #take the longitude and latitude data and clean it
  gis_data = family_location_data[3:4]
  gis_data = na.omit(gis_data)
  gis_data = unique(gis_data)

  #define bound region of raster map to speed up data processing
  lon_min = min(gis_data$longitude)
  lon_max = max(gis_data$longitude)
  lat_min = min(gis_data$latitude)
  lat_max = max(gis_data$latitude)

  #get the requested BioClim data
  bio_clim_data = raster(paste0("bio_", zone_number, ".bil"))
  #check to see if the BioClim data has been correctly extracted
  #fromDisk(bio_clim_data)

  #crop BioClim map to make it easier to precess
  bio_clim_data = crop(bio_clim_data, extent(c(lon_min, lon_max, lat_min, lat_max)))

  #correct data: it is stored as a decimal ten times larger than actual to maintain
  #precision and to save memory
  bio_clim_data <- bio_clim_data/10

  #make the position data have the correct longitude and latitude format for
  #the projection used by the raster data
  #str(gis_data)
  coordinates(gis_data) <- ~longitude+latitude
  #"The data is in the latitude/longitude coordinate reference system (not projected) and the datum is WGS84"
  #proj4string(gis_data) <- CRS("+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs")

  #extract BioClim data at each reported locations
  bio9dat=extract(bio_clim_data, gis_data)

  #if user want to store in the entered directory
  if (storeInDirectory) {
    #output
    if (!is.null(filename)) {
      write.csv(bio9dat, file = filename)
    }
    else {
      write.csv(bio9dat, file=paste0("bio_", zone_number, "_data.csv"))
    }
  }
  #otherwise store in old directory
  else {
    setwd(old_directory)
    #output
    if (!is.null(filename)) {
      write.csv(bio9dat, file = filename)
    }
    else {
      write.csv(bio9dat, file=paste0("bio_", zone_number, "_data.csv"))
    }
  }

  #restore working directory to original directory regardless of save preference and return data
  setwd(old_directory)
  return(bio9dat)
}






#' Raster Layer Coordinates
#'
#' @param longitude A number -180 to 180
#' @param latitude A number -90 to 90
#' @param raster_layer A raster layer
#'
#' @export
#'
#' @return This function returns the values of a raster layer at the specified coordinates.
RasterLayerCoordinates <- function(longitude, latitude, raster_layer) {
  #invoke required packages
  require(sp)
  require(raster)

  #if there is an unequal number of longitude and latitude points
  #something is wrong on the function call end and so we stop function
  try(if(length(longitude)!=length(latitude)) stop("Unequal number of values in longitude and latitude."))

  #build data frame
  coordinate_pairs = data.frame(lon = longitude, lat = latitude)

  #convert data frame to spacial data frame
  coordinates(coordinate_pairs) <- ~lon+lat

  #extract data
  extracted_values = extract(raster_layer, coordinate_pairs)

  return (extracted_values)
}
