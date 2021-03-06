% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ExtractingMultipleRasterLayers.R
\name{ExtractMultiple}
\alias{ExtractMultiple}
\title{Extracting from Multiple Raster Layers}
\usage{
ExtractMultiple(raster_vector, gis_data, projection = NULL, filename = NULL,
  storeInDirectory = FALSE)
}
\arguments{
\item{raster_vector}{A single raster layer or a vector}

\item{gis_data}{A data frame containing longitude and latitude pairs}

\item{projection}{User specified projection #may not do anything becasue extract may reset it}

\item{filename}{A user specified name that the data will be saved as, instead of the default raster_data.csv}

\item{storeInDirectory}{This boolean parameter asks wheter or not the user wants to store
the data in the working directory.}
}
\value{
This function returns a data frame/matrix at each location specified by gis_data parameter.
If storeInDirectory is TRUE, then it also saves the data to a file (defaults to raster_data.csv).
}
\note{
I'm using a raster brick because, "while a RasterBrick has to refer to one multi-layer file or
is in itself a multi-layer object with data loaded in memory, a RasterStack may 'virtually'
connect several raster objects written to different files or in memory. Processing will be
more efficient for a RasterBrick than for a RasterStack, but RasterStack has the advantage
of facilitating pixel based calculations on separate raster layers."
From http://geoscripting-wur.github.io/IntroToRaster/
}
\examples{
ExtractMultiple(c(bio1, bio9), gis_data)
species_vector<-c("Abies amabilis", "Acer nigrum")
gis_data = BIEN.gis.species(species_vector)[2:3]
gis_data = na.omit(gis_data)
gis_data = unique(gis_data)
raster_vec <- c(raster("bio_1.bil"), raster("bio_9.bil"), raster("bio_4.bil"))
ExtractMultiple(raster_vector = raster_vec, gis_data = gis_data, filename = "your_file_name.csv", storeInDirectory = TRUE)
}
\author{
Alton Barbehenn
}
\keyword{extract}
\keyword{layers}
\keyword{raster}

