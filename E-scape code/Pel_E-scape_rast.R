#' """ Make E-scape for menhaden for Pelicans
#'     Uses random points to generate IEI
#'     Made from cover raster and chla raster
#'     @author = Ryan James
#'     Date = 6/23/21"""

#setwd("C:/Users/brudi/Nelson Lab Dropbox/James Nelson/MenEscape")

setwd("C:/Users/Wjames/OneDrive - Florida International University/PelE-scape")

library(raster)
library(landscapemetrics)
library(tidyverse)
library(sf)
library(stars)
library(exactextractr) 
library(fasterize)
library(terra)

setwd("C:/Users/Wjames/OneDrive - Florida International University/PelE-scape")

r = rast('gis/LAmarsh13.tif') 

# scale raster for 17-19
chl = rast('gis/chlaMean17-19.tif')
crs(chl) = "epsg:4326"
chl = project(chl, crs(r))
setMinMax(chl)
chl = chl/minmax(chl)[2]
writeRaster(chl, 'gis/chlaScale17-19.tif')

# scale raster for 19-21
chla = rast('gis/chlaMean19-21.tif')
crs(chla) = "epsg:4326"
chla = project(chla, crs(r))
setMinMax(chla)
chla = chla/minmax(chla)[2]
writeRaster(chla, 'gis/chlaScale19-21.tif')

# extend raster 
r = rast('gis/LAmarsh13.tif') 
e = ext(399112.499999412, 950002.49999941, 3050000, 3385268.49990729)
re = extend(r, e, fill = 6)
writeRaster(re, 'gis/LAmarsh13E.tif')
