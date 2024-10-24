#' """ menhaden MixSIAR w/ marine and marsh only
#'    @author: Ryan James
#'    Date: 1/25/18
#'    Editted: 1/19/21"""

setwd("C:/Users/wrjam/Dropbox/WorkDocs/R/MenEscape")
setwd("C:/Users/Ryan/Dropbox/WorkDocs/R/MenEscape")
library(MixSIAR)
options(max.print=6000000)

#### Menhaden with marine POM ****this is the one****----
## load data
mix = load_mix_data(file("data/Menhaden.csv"),
                     iso_names=c("d13C","d34S"),
                     factors= c('Species'),
                     fac_random=c(F),
                     fac_nested=c(FALSE),
                     cont_effects=NULL)

source = load_source_data(file("data/sMenZoCorrected.csv"),
                           source_factors=NULL,
                           conc_dep=FALSE,
                           data_type="means",
                           mix)

discr = load_discr_data(file("data/TEF.csv"), mix)

#  # Make an isospace plot
plot_data(filename="isospace_plot", plot_save_pdf=FALSE, plot_save_png=FALSE, mix,source,discr)

# Write the JAGS model file
model_filename = "MixSIAR_model.txt"
resid_err = T
process_err = T
write_JAGS_model(model_filename, resid_err, process_err, mix, source)

jags.men = run_model(run="very long", mix, source, discr, model_filename,
                      alpha.prior = 1, resid_err, process_err)

# Process JAGS output
output_men = list(summary_save = TRUE,
                   summary_name = "data/menCS_ss",
                   sup_post = FALSE,
                   plot_post_save_pdf = FALSE,
                   plot_post_name = "mid_posterior_density",
                   sup_pairs = FALSE,
                   plot_pairs_save_pdf = FALSE,
                   plot_pairs_name = "mid_pairs_plot",
                   sup_xy = TRUE,
                   plot_xy_save_pdf = FALSE,
                   plot_xy_name = "mid_xy_plot",
                   gelman = TRUE,
                   heidel = FALSE,
                   geweke = TRUE,
                   diag_save = T,
                   diag_name = "data/menCS_diag",
                   indiv_effect = FALSE,
                   plot_post_save_png = FALSE,
                   plot_pairs_save_png = FALSE,
                   plot_xy_save_png = FALSE)

output_JAGS(jags.men, mix, source, output_men)



# run function to output csv of MixSIAR results
ind = mixTable('menmarine_ss.txt', 'With Marine POM')
write.csv(ind,'menmixmarine.csv', row.names = F)

#### Menhaden with out marine POM----
## load data
mix = load_mix_data(file("data/Menhaden.csv"),
                    iso_names=c("d13C","d15N","d34S"),
                    factors= c('Species'),
                    fac_random=c(F),
                    fac_nested=c(F),
                    cont_effects=NULL)

source = load_source_data(file("data/sMenMarineConc.csv"),
                          source_factors=NULL,
                          conc_dep=FALSE,
                          data_type="means",
                          mix)

discr = load_discr_data(file("data/TEFmarsh2.csv"), mix)

#  # Make an isospace plot
plot_data(filename="isospace_plot", plot_save_pdf=FALSE, plot_save_png=FALSE, mix,source,discr)

# Write the JAGS model file
model_filename = "MixSIAR_model.txt"
resid_err = T
process_err = T
write_JAGS_model(model_filename, resid_err, process_err, mix, source)

jags.menmarsh = run_model(run="normal", mix, source, discr, model_filename,
                     alpha.prior = 1, resid_err, process_err)

# Process JAGS output
output_menmarsh = list(summary_save = TRUE,
                  summary_name = "data/menmarsh_ss",
                  sup_post = FALSE,
                  plot_post_save_pdf = FALSE,
                  plot_post_name = "mid_posterior_density",
                  sup_pairs = FALSE,
                  plot_pairs_save_pdf = FALSE,
                  plot_pairs_name = "mid_pairs_plot",
                  sup_xy = TRUE,
                  plot_xy_save_pdf = FALSE,
                  plot_xy_name = "mid_xy_plot",
                  gelman = TRUE,
                  heidel = FALSE,
                  geweke = TRUE,
                  diag_save = FALSE,
                  diag_name = "mid_diagnostics",
                  indiv_effect = FALSE,
                  plot_post_save_png = FALSE,
                  plot_pairs_save_png = FALSE,
                  plot_xy_save_png = FALSE)

output_JAGS(jags.menmarsh, mix, source, output_menmarsh)
# run function to output csv of MixSIAR results
MenMarsh = mixTable('menmarsh_ss.txt', 'Marsh only')
write.csv(MenMarsh,'menmixmarsh.csv', row.names = F)


#### Menhaden Average of sp with marine POM----
## load data
mix = load_mix_data(file("data/Menhaden.csv"),
                    iso_names=c("d13C","d15N","d34S"),
                    factors= c('Species'),
                    fac_random= F,
                    fac_nested= F,
                    cont_effects=NULL)

source = load_source_data(file("data/sMenhadenConc.csv"),
                          source_factors=NULL,
                          conc_dep=FALSE,
                          data_type="means",
                          mix)

discr = load_discr_data(file("data/TEFpom.csv"), mix)

#  # Make an isospace plot
plot_data(filename="isospace_plot", plot_save_pdf=FALSE, plot_save_png=FALSE, mix,source,discr)

# Write the JAGS model file
model_filename = "MixSIAR_model.txt"
resid_err = T
process_err = T
write_JAGS_model(model_filename, resid_err, process_err, mix, source)

jags.men = run_model(run="very long", mix, source, discr, model_filename,
                     alpha.prior = 1, resid_err, process_err)

# Process JAGS output
output_men = list(summary_save = TRUE,
                  summary_name = "data/menone_ss",
                  sup_post = FALSE,
                  plot_post_save_pdf = FALSE,
                  plot_post_name = "mid_posterior_density",
                  sup_pairs = FALSE,
                  plot_pairs_save_pdf = FALSE,
                  plot_pairs_name = "mid_pairs_plot",
                  sup_xy = TRUE,
                  plot_xy_save_pdf = FALSE,
                  plot_xy_name = "mid_xy_plot",
                  gelman = TRUE,
                  heidel = FALSE,
                  geweke = TRUE,
                  diag_save = FALSE,
                  diag_name = "mid_diagnostics",
                  indiv_effect = FALSE,
                  plot_post_save_png = FALSE,
                  plot_pairs_save_png = FALSE,
                  plot_xy_save_png = FALSE)

output_JAGS(jags.men, mix, source, output_men)

mixTable = function(file,site){
  require(tidyverse)
  b = read_file(file)
  z = (str_split_fixed(b,pattern = "M", n=2)[1] %>% str_length()) -
    (str_split_fixed(b,pattern = "M", n=2)[1] %>% str_trim('right') %>% str_length())- 5
  s = c(z,6,6,6,6,6,6,6,6,6)
  cn = c('ID', 'Mean', 'SD', '2.5%', '5%', '25%', '50%', '75%', '95%', '97.5%')
  x = read_fwf(file, skip = 8, fwf_widths(s, cn))
  x$source = 0
  x$name = 0
  for (i in 1:nrow(x)){
    temp = strsplit(x$ID, split = '.', fixed = T)
    x$source[i] = temp[[i]][3]
    x$name[i] = temp[[i]][2]
  }
  
  x$site = site
  
  df = data.frame(x$name, x$site, x$source, x$Mean, x$SD, x$`2.5%`, x$`97.5%`)
  colnames(df) = c('name', 'site', 'source', 'mean', 'sd', 'lowend', 'highend')
  df = na.omit(df)
  
  return(df)
}

# run function to output csv of MixSIAR results
ind = mixTable('menone_ss.txt', 'average POM')
write.csv(ind,'menmixone.csv', row.names = F)


# figures for mm -----

# function to convert MixSIAR output to a table to generate random data
# to use for hypervolumes
# file = name of .txt of summary statistics of MixSIAR
# site = 'site' or 'location' of the data set
mixTable = function(file,site,ind = NULL, global = T){
  require(tidyverse)
  cn = c('ID', 'Mean', 'SD', '2.5%', '5%', '25%', '50%', '75%', '95%', '97.5%')
  x = read_table(file, skip = 7,  col_names = T)
  names(x) = cn
  x$source = NA
  x$name = NA
  
  for (i in 1:nrow(x)){
    temp = strsplit(x$ID, split = '.', fixed = T)
    x$source[i] = temp[[i]][3]
    x$name[i] = temp[[i]][2]
  }
  
  x$site = site
  x$ymax = x$`75%` + 1.5*(x$`75%` - x$`25%`)
  x$ymin = x$`25%` - 1.5*(x$`75%` - x$`25%`)
  
  df = data.frame(x$name, x$site, x$source, x$Mean, x$SD, x$`2.5%`, x$`97.5%`,
                  x$`50%`, x$`25%`, x$`75%`, x$ymax, x$ymin)
  colnames(df) = c('name', 'site', 'source', 'mean', 'sd', 'lowend', 'highend',
                   'mid', 'low', 'up', 'ymax', 'ymin')
  for (i in 1:nrow(df)){
    if (df$ymax[i] > df$highend[i]){
      df$ymax[i] = df$highend[i]
    }
    if (df$ymin[i] < df$lowend[i]){
      df$ymin[i] = df$lowend[i]
    }
  }
  df = na.omit(df)
  if(global == T){
   df = subset(df, name != 'global') 
   }
  if (isTRUE(ind == T) == T){
    df = data.frame(x$name, x$site, x$source, x$Mean)
    colnames(df) = c('name', 'site', 'source', 'mean')
    df = spread(df, 'source', 'mean')
  }
  
  return(df)
}

d = mixTable('menone_ss.txt', 'multiple POM', global = F)

db = rbind(d[1,], d[4,])
db = rbind(db, b)
ggplot(db, aes(x = source,  fill = source))+
  geom_boxplot(aes(middle = mid,
                   upper = up,
                   lower = low,
                   ymax = ymax,
                   ymin = ymin), stat = 'identity')+
  theme_classic()+ 
  scale_x_discrete(labels = c('BMA', expression(italic('Spartina')), 'Phytoplankton'))+
  scale_fill_manual(values = c('tan4', 'darkgreen', 'deepskyblue3'))+
  theme(legend.position = 'none')+
  labs(x = NULL, y = 'Source contribution')+
  theme(axis.title = element_text(size = 20), 
        axis.text = element_text(size = 20, colour = "gray0"), 
        plot.title = element_text(size = 20))
ggsave("menmm.tiff", units="in", width=8, height=6, dpi=600,compression = 'lzw') 


