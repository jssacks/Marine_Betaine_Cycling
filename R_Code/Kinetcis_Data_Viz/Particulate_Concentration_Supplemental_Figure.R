

#packages:
library(tidyverse)
library(ggpubr)



#define inputs:
flux.kin.file <- "Intermediates/Compiled_Kinetics_Fluxes.csv"
context.file <- "Intermediates/flux_contextualizaing_data"
meta.data.file <- "Intermediates/Kinetics_Meta_Data_Compiled.csv"
enviro.metab.file <- "Intermediates/Kinetics_Exp_Enviro_Metabolome_Dat.csv"

region.pal <- c("#00887d", "#014d64", "#01a2d9", "#6794a7")



#load in and combine datasets:

#load in data:

#region classifications
region.dat <- read_csv(meta.data.file) %>%
  select(KinExp_ID, Region) %>%
  separate(KinExp_ID, into = c("Cruise", "exp"), remove = FALSE) %>%
  mutate(Region = case_when(Region == "Gyre" ~ "NPSG",
                            Region == "Coastal" ~ "Puget Sound",
                            TRUE ~ Region)) %>%
  mutate(Region = as.factor(Region)) %>%
  mutate(Region = fct_relevel(Region, c("NPTZ", "NPSG", "Equatorial", "Puget Sound"))) 

##Remove TMAO and calculate total betaine + sulfonium concentration:
summed.comp.dat <- read_csv(enviro.metab.file) %>%
  filter(!Compound == "TMAO") %>%
  group_by(Cruise, exp) %>%
  reframe(metab.tot.nM = sum(Mean.Diss.Conc.nM),
          metab.tot.nM.uncert = sqrt(sum(SD.Diss.Conc.nM)^2))

#flux and kinetics data
flux.kin.dat <- read_csv(flux.kin.file) %>%
  select(Cruise, exp, Compound, Mean.Diss.Conc.nM, SD.Diss.Conc.nM, mean_ks, sd_ks, mm_flux_nM_day, mm_flux_sd_nM_day)

#particulate concentration data
context.dat <- read_csv(context.file) %>%
  select(KinExp_ID, Cruise, exp, Compound, mean.part.conc.nM, sd.part.conc.nM)


##Combine all datasets
all.dat <- flux.kin.dat %>%
  left_join(., context.dat) %>%
  left_join(., summed.comp.dat) %>%
  left_join(., region.dat) 


###Make particulate concentration figure:
p.plot <- ggplot(all.dat, aes(x = KinExp_ID, y = mean.part.conc.nM, fill = Region)) +
  geom_errorbar(aes(ymin = mean.part.conc.nM-sd.part.conc.nM, ymax = mean.part.conc.nM+sd.part.conc.nM), width = 0.2) +
  geom_point(size = 3, shape = 21) +
  facet_grid(.~Compound, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = region.pal) +
  scale_y_log10() +
  theme_bw() +
  xlab("Experiment") +
  ylab("Particulate Concentration (nM)") +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.8, hjust = 0.8),
        #  strip.background.x = element_blank(),
        strip.text.x = element_text(face = "bold")) 
p.plot

ggsave(p.plot, filename = "Figures/Flux_Paper_Figures/Particulate_Metabolite_Concentration_Plot.png",
       dpi = 600, units = "in", height = 3.5, width = 5, scale = 1.2)
