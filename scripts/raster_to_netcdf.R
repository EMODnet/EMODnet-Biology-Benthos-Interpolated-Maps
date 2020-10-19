# Run in ../data as:
# Rscript ../script/raster_to_netcdf.R

# get data/EUSM2019_substrate.gdb from EMODnet Seabed habitat

# install.package("sf")
# install.packages("raster")
# install.packages("fasterize")
# install.packages("ncdf4")
# install.packages("RNetCDF")

require(sf)
require(raster)
require(fasterize)

library(ncdf4)
library(RNetCDF)

seds<-sf::st_read(dsn = "EUSM2019_substrate.gdb", layer = "EUSM2019_Arctic_Atlantic_substrate")
#r<-raster(seds,res=10000)
#bb<-fasterize(seds,r,by="Substrate")

seds_lonlat <- st_transform(seds, "+proj=longlat +datum=WGS84 +no_defs")
r<-raster(seds_lonlat,res=0.05)
bb<-fasterize(seds_lonlat,r,by="Substrate")

writeRaster(bb, "substrate.nc", overwrite=TRUE, format="CDF", bylayer=TRUE,suffix='names')

