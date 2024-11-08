#' """ E-scape plot for menhaden
#'     @author = Ryan James
#'     Date = 7/9/23"""

library(tidyverse)
library(sf)
library(stars)
library(ggpubr)

setwd("C:/Users/wrjam/OneDrive - Florida International University/PelE-scape")

m19 = read_stars('E-scapes/menE-scape17-19_1km.tiff') |> 
  st_as_sf() |> 
  rename(HRI = `menE-scape17-19_1km.tiff`)

m21 = read_stars('E-scapes/menE-scape19-21_1kn.tiff') |> 
  st_as_sf() |> 
  rename(HRI = `menE-scape19-21_1kn.tiff`)

la = st_read('gis/Gulf_coast_states.shp') 

a = ggplot() + 
  geom_sf(data = la, fill = '#DEDEDE', color = NA)+
  geom_sf(data = m19, aes(fill = HRI), color = NA)+
  theme_bw()+
  labs(title = expression(italic(E)*'-scape 2019'))+
  scale_fill_gradientn(colors = c("#0571b0",'floralwhite',"#b5001c"),
                       limits = c(0,7.1),
                       values=scales::rescale(c(0,1,7.1)),
                       breaks = c(0,1,2,3,4,5,6,7))+
  coord_sf(ylim = c(28.7, 30.5), xlim = c(-93.8,-88.5), expand = F)+
  theme(axis.title = element_text(size = 18), 
        axis.text = element_text(size = 14, colour = "gray0"), 
        #panel.background = element_rect(fill="aliceblue"),
        plot.title = element_text(size = 18, hjust=0.5),
        legend.title = element_text(size = 14),
        strip.text.x = element_text(size = 18),
        legend.text = element_text(size = 13))

b = ggplot() + 
  geom_sf(data = la, fill = '#DEDEDE', color = NA)+
  geom_sf(data = m21, aes(fill = HRI), color = NA)+
  theme_bw()+
  labs(title = expression(italic(E)*'-scape 2021'))+
  scale_fill_gradientn( colors = c("#0571b0",'floralwhite',"#b5001c"),
                        limits = c(0,7.1),
                        values=scales::rescale(c(0,1,7.1)),
                        breaks = c(0,1,2,3,4,5,6,7))+
  coord_sf(ylim = c(28.7, 30.5), xlim = c(-93.7,-88.5), expand = F)+
  theme(axis.title = element_text(size = 18), 
        axis.text = element_text(size = 14, colour = "gray0"), 
        #panel.background = element_rect(fill="aliceblue"),
        plot.title = element_text(size = 18, hjust=0.5),
        legend.title = element_text(size = 14),
        strip.text.x = element_text(size = 18),
        legend.text = element_text(size = 13))

aos.escapes <- ggarrange(a, b, ncol = 1)

ggsave('figs/E-scapeA_aos.tif', plot = a, dpi = 300, units = 'in', 
       width = 9, height = 4.5)

ggsave('figs/E-scapeB_aos.tif', plot = b, dpi = 300, units = 'in', 
       width = 9, height = 4.5)
