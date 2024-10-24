#' """ Make E-scape for menhaden for Pelicans
#'     Uses random points to generate IEI
#'     Made from cover raster and chla raster
#'     @author = Ryan James
#'     Date = 6/23/21
#'     Edit = 1/8/23"""

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


# function to calculate area and edge habitat around sampling point
# r = base raster file of sampling area 
# chla = chlorophyll a raster scaled so max value = 1
# data = tibble with mixing model results for a single species ***has to be one row***
# radius = radius around point to calculate metrics
# n = number of points to randomly sample
randIEI = function(r, chla, data, radius, n){
  require(raster)
  require(exactextractr)
  require(landscapemetrics)
  require(tidyverse)
  require(sf)
  require(terra)
  
  # make polygon of smaller raster 
  poly = st_as_sfc(st_bbox(r))
  
  # random points until n samples met on raster with values
  samp = st_sample(poly, n*50) %>%
    st_as_sf 
  
  # extract values of chla for samp points
  ext = samp %>% 
    mutate(ext = terra::extract(chla, samp)) %>%
    drop_na
  
  # grab n point at random from points with values
  samplocs = sample_n(ext, n, replace = FALSE)%>%
    select(- ext)
  
  # bind mixing model information
  buf = cbind(samplocs, data) %>%
    st_transform(st_crs(r)) %>%# utm 15
    mutate(site = row_number())%>%
    st_buffer(dist = radius)
  
  # make empty data frame to calculate cover area at each site
  df = tibble(site = buf$site, f_marsh = NA, f_water = NA, 
              f_edge = NA, f_chla = NA, f_waterChla = NA)
  
  # calculate the cover areas around sampling point in circle of size radius
  # in the habitat 
  met = sample_lsm(r, y = buf, plot_id = buf$site, 
                   what = c('lsm_c_pland'), 
                   return_raster = T, progress = F)
                    
  
  # assign cover area and measure edge area
  for (i in 1:nrow(df)){
    
    # filter out single plot
    a = met %>% filter(plot_id == df$site[i], metric == 'pland') 
    
    # initilize the areas of habitat types
    f_marsh = 0
    f_water = 0
    f_edge = 0 
    
    if (a$percentage_inside[1] > 0) {
    # renames cover types from 1 and 6 to marsh, water
   
      for (j in 1:nrow(a)){
        if (a$class[j] == 1){
         f_marsh = a$value[j]/100
        }else if (a$class[j] == 6){
         f_water = a$value[j]/100
        }
      } 
    
    
     # calculate the total edge between habitat classes
     # multiply the number of cells that are adjacent
     # and multiply by resolution of cells
     per = get_adjacencies(a$raster_sample_plots[1])[[1]]*res(a$raster_sample_plots[[1]])[1]
    
     # calculate total area of raster 
     t = lsm_l_ta(a$raster_sample_plots[1])
     ta = t$value*10000
    
     # calculate the edge of marsh (class 1) and water (class 6)
     # multiply
     if (f_water > 0 & f_marsh > 0){
       mar_edge = per['1','6']*res(a$raster_sample_plots[[1]])[1]
     }else{
       mar_edge = 0
     }
    
     # calculate f_edge
     f_edge = (mar_edge)/ta
    }
    
    df$f_marsh[i] = f_marsh
    df$f_water[i] = f_water
    df$f_edge[i] = f_edge
    
  }
  # assign chla mean for each buffer
  df = df %>% 
    mutate(f_chla = exact_extract(chla, buf, 'mean'))
  
  # calculate water production
  df = df %>% 
    mutate(f_waterChla = f_water * f_chla)
  
  # bind random points together
  cd = st_transform(samplocs, 4326)
  ll = unlist(st_geometry(cd)) %>% 
    matrix(ncol=2,byrow=TRUE) %>% 
    as_tibble() %>% 
    setNames(c("lon","lat"))
  
  pts = buf%>% as_tibble %>% cbind(ll) %>% select(-x)
  
  # calculate IEI assigning NA if none of that habitat type
  d = full_join(pts, df, by = 'site')%>% as_tibble() %>% 
    mutate(edgeiei = BA/na_if(f_edge,0),
           marshiei = SPART/na_if(f_marsh,0),
           wateriei = POM/na_if(f_waterChla,0))

  
  return(d)
}

# function to make an E-scape
# r = raster of area
# iei =tibble output from randIEI function 
# radius = radius of circle used to calculate IEI and *2 will be length of cell to generate E-scape
# shp = shape file to clip the E-scape polygon grid
# rast = T or F; if true return raster of E-scape. If false returns sf
# prog = T or F; default false, if true show progress of landscape metrics calcs
E_scape = function(r, chla, iei, radius, shp = NA, rast = F, prog = F){
  require(raster)
  require(landscapemetrics)
  require(tidyverse)
  require(sf)
  require(fasterize)
  require(exactextractr)
  
  # make grid over raster
  e = ext(r)
  bb = st_bbox(c(e[1],
                 e[3],
                 e[2],
                 e[4]), crs = st_crs(r)) 
  grid_geom = st_make_grid(st_as_sfc(bb), cellsize = radius*2)
  grid = st_sf(geom = grid_geom) %>% 
    mutate(site = row_number())
  
  if(!is.na(shp)){
    grid =  grid %>% filter(lengths(st_intersects(., shp)) > 0)
  }
  # calculate the cover areas within grid
  met = sample_lsm(r, grid, plot_id = grid$site,
                   level = "class", metric = 'pland',
                   return_raster = T, progress = prog)
  
  # make empty data frame to calculate cover area at each site
  df = tibble(site = grid$site, f_marsh = NA, f_water = NA, 
              f_edge = NA, f_chla = NA, f_waterChla = NA, pi = NA)
  
  # assign cover area and measure edge area
  for (i in 1:nrow(df)){
    
    # filter out single plot
    a = met %>% filter(plot_id == df$site[i], metric == 'pland') 
    
    # initilize the areas of habitat types
    f_marsh = 0
    f_water = 0
    f_edge = 0 
    
    if (a$percentage_inside[1] > 0) {
      # renames cover types from 1 and 6 to marsh, water
      
      for (j in 1:nrow(a)){
        if (a$class[j] == 1){
          f_marsh = a$value[j]/100
        }else if (a$class[j] == 6){
          f_water = a$value[j]/100
        }
      } 
      
      
      # calculate the total edge between habitat classes
      # multiply the number of cells that are adjacent
      # and multiply by resolution of cells
      per = get_adjacencies(a$raster_sample_plots[1])[[1]]*res(a$raster_sample_plots[[1]])[1]
      
      # calculate total area of raster 
      t = lsm_l_ta(a$raster_sample_plots[1])
      ta = t$value*10000
      
      # calculate the edge of marsh (class 1) and water (class 6)
      # multiply
      if (f_water > 0 & f_marsh > 0){
        mar_edge = per['1','6']*res(a$raster_sample_plots[[1]])[1]
      }else{
        mar_edge = 0
      }
      
      # calculate f_edge
      f_edge = (mar_edge)/ta
    }
    
    df$f_marsh[i] = f_marsh
    df$f_water[i] = f_water
    df$f_edge[i] = f_edge
    df$pi[i] = a$percentage_inside[1]
    
    if (i/nrow(df) == 0.05){
      cat('5% done \n')
    } else if (i/nrow(df) == 0.25) {
      cat('25% done \n')
    } else if (i/nrow(df) == 0.5) {
      cat('50% done \n')
    } else if (i/nrow(df) == 0.75){
      cat('75% done \n')
    } else if (i/nrow(df) == 0.95){
      cat('95% done \n')
    }
    
    
  }
  # assign chla mean for each buffer
  df = df %>% 
    mutate(f_chla = exact_extract(chla, grid, 'mean'))
  
  # calculate water production
  df = df %>% 
    mutate(f_waterChla = f_water * f_chla)
  
  # calculate HRI
  edgeiei = median(iei$edgeiei, na.rm = T)
  marshiei = median(iei$marshiei, na.rm = T)
  wateriei = median(iei$wateriei, na.rm = T)
  
  df = df %>% mutate(
    HRI = f_marsh*marshiei + f_edge*edgeiei + f_waterChla*wateriei
  )
 
  d = full_join(grid, df, by = 'site')
  
  if(!is.na(shp)){
    d = st_intersection(d, shp)
    d = st_collection_extract(d, 'POLYGON')
    #d = d %>% group_by(Site) %>% summarise_all(mean)
    
  }
  
  
  if (rast == T){
    x = (extent(d)[2]-extent(d)[1]) %/% (radius*2)
    y = (extent(d)[4]-extent(d)[3]) %/% (radius*2)
    
    # make a raster of the points from the area
    ras = raster(ncol = x, nrow = y)
    extent(ras) = extent(d)
    crs(ras) = crs(r)
    
    
    d = fasterize(d, ras, field = 'HRI')
  }
  
  return(d)
} 

# Calculate IEI ----
# load rasters
#r = rast('gis/LAmarsh13E.tif') 
r = rast('gis/LAmarsh13.tif') 
chl = rast('gis/chlaScale17-19.tif')

chla = raster('gis/chlaScale19-21.tif')

# mixing model results
data = tibble(BA = 0.228, POM = 0.680, SPART = 0.092)

# size of foraging in m
radius = 500

# number of random points 
n = 1000

# rand IEI
# 17-19
menIEI = randIEI(r, chl, data, radius, n)

write_csv(menIEI, 'data/BGmenIEI17-19_n1000_r500.csv')


# 19-21
menIEI = randIEI(r, chla, data, radius, n)
write_csv(menIEI, 'data/BGmenIEI19-21_n1000_r500.csv')

# marshIEI = median(menIEI$marshiei, na.rm = T)
# edgeIEI = median(menIEI$edgeiei, na.rm = T)
# waterIEI = median(menIEI$wateriei, na.rm = T)

#17-19
iei = read_csv('data/BGmenIEI17-19_n1000_r500.csv')

menE_scape = E_scape(r, chl, iei, radius, prog = F)

m = menE_scape %>% mutate(HRI = na_if(HRI, 0))

write_csv(m, 'data/BGmenE-scape17-19_500m.csv')

# convert SF objects into rasters
# 17-19

# convert to raster
x = (extent(menE_scape)[2]-extent(menE_scape)[1]) %/% (radius*2)
y = (extent(menE_scape)[4]-extent(menE_scape)[3]) %/% (radius*2)

# make a raster of the points from the area
ras = raster(ncol = x, nrow = y)
extent(ras) = extent(menE_scape)
crs(ras) = crs(r)


d = fasterize(m, ras, field = 'HRI')
writeRaster(d, 'E-scapes/BGmenE-scape17-19_1km.tiff', overwrite = T)

#19-21
iei = read_csv('data/BGmenIEI19-21_n1000_r500.csv')

menE_scape = E_scape(r, chla, iei, radius, prog = F)

m = menE_scape %>% mutate(HRI = na_if(HRI, 0))

write_csv(m, 'data/BGmenE-scape19-21_500m.csv')


# convert SF objects into rasters
# 17-19

# convert to raster
x = (extent(menE_scape)[2]-extent(menE_scape)[1]) %/% (radius*2)
y = (extent(menE_scape)[4]-extent(menE_scape)[3]) %/% (radius*2)

# make a raster of the points from the area
ras = raster(ncol = x, nrow = y)
extent(ras) = extent(menE_scape)
crs(ras) = crs(r)


d = fasterize(m, ras, field = 'HRI')
writeRaster(d, 'E-scapes/BGmenE-scape19-21_1kn.tiff', overwrite = T)

ggplot() + 
  geom_sf(data = m, aes(fill = log(HRI)), color = NA)+
  theme_bw()+
  scale_fill_gradientn( colors = c("#0571b0",'floralwhite',"#ca0020"),
                        values=scales::rescale(c(-9,0,2)))+
  theme(axis.title = element_text(size = 18), 
        axis.text = element_text(size = 14, colour = "gray0"), 
        #panel.background = element_rect(fill="aliceblue"),
        plot.title = element_text(size = 18, hjust=0.5),
        legend.title = element_text(size = 14),
        strip.text.x = element_text(size = 18),
        legend.text = element_text(size = 13))




read_csv('data/BGmenIEI17-19_n1000_r500.csv') |> 
  reframe(across(edgeiei:wateriei, \(x) quantile(x, c(0.025, 0.25, 0.5, 0.75, 0.975), 
                                                 na.rm = T))) 

read_csv('data/BGmenIEI19-21_n1000_r500.csv') |> 
  reframe(across(edgeiei:wateriei, \(x) quantile(x, c(0.025, 0.25, 0.5, 0.75, 0.975), 
                                                 na.rm = T))) 


read_csv('data/BGmenIEI19-21_n1000_r500.csv') |> 
  bind_rows(read_csv('data/BGmenIEI17-19_n1000_r500.csv')) |> 
  reframe(across(edgeiei:wateriei, \(x) quantile(x, c(0.025, 0.25, 0.5, 0.75, 0.975), 
                                                 na.rm = T))) 


read_csv('data/BGmenIEI17-19_n1000_r500.csv') |> 
  reframe(across(edgeiei:wateriei, \(x) quantile(x, c(0.5, 0.25, 0.75), 
                                                 na.rm = T))) 

read_csv('data/BGmenIEI19-21_n1000_r500.csv') |> 
  reframe(across(edgeiei:wateriei, \(x) quantile(x, c(0.5, 0.25, 0.75), 
                                                 na.rm = T))) 
