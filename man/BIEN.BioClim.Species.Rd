% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BioClimData.R
\name{BIEN.BioClim.Species}
\alias{BIEN.BioClim.Species}
\title{A Bio Clim data for a vector of species}
\usage{
BIEN.BioClim.Species(species_vector, zone_number, filename = NULL,
  directory = NULL, storeInDirectory = FALSE)
}
\arguments{
\item{species_vector}{A vector of species in the BIEN3 database. This should be a vector of strings, for example species_vector<-c("Abies amabilis", "Acer nigrum").}

\item{zone_number}{The zone number is an integer 1-18 (inclusive) representing that specific Bio Clim, for example, Bio Clim 4 would be zone_number = 4.}

\item{filename}{The location of the exported data, if you choose to use that parameter. Defaults to NULL.}

\item{directory}{Location of the Bio Clim data. Defaults to NULL.}

\item{storeInDirectory}{A boolean parameter representing whether or not the exported data should be in the working directory (FALSE) or with the Bio Clim data (TRUE). Defaults to FALSE.}
}
\value{
This function returns the value(s) of a given Bio Clim raster layer at every reported point for a species vector.
This function exports a csv containing the same information that is returned.
}
\description{
This function is utilizes functions from the BIEN3 Project Written by Brian Maitner 2015.
It takes as input a vector containing a set of species and a zone number. This function
then queries the BIEN3 database for all reported locations of a memeber of the species.
Then it creates a csv file containing the bioclim values of all these locations.
}
\details{
As an optional parameters, the user can choose to enter thier own file name for the csv file,
or let the code use the default name. Also if you wish to use a .bil file in a different
the user can enter this as a string. Bio Clim data must be in the working directory or in the
directoery specified by the directory parameter.
}
\note{
The examples are labeled as "Not run" becuase the function accesses files and if the directory
is not correct, when the code is tested it wont be able to find the .bil file.
}
\examples{
\dontrun{BIEN.BioClim.Species(species_vector = s_vec, zone_number = 7)
species_data <- BIEN.BioClim.Species(species_vector = s_vec, zone_number = 7, filename = "test.csv", directory = "/filelocation/", storeInDirectory = TRUE)}

}

