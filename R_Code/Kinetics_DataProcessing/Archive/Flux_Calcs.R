

####
library(tidyverse)






#####Define inputs:
tt.file <- "Intermediates/turnover_time_estimates.csv"
dissolved.file <- "Intermediates/Dissolved_Quantified_Data.csv"
match.file <- "Meta_Data/Ingalls_Lab_Data/KinExp_DissSamp_Match_Key.csv"

  
## tt dat
tt.dat <- read_csv(tt.file) %>%
  rename("KinExp_ID" = exp) 

#match file
match.dat <- read_csv(match.file) %>%
  mutate(Compound = str_replace(Compound, "GBt", "Glycine betaine"))
repid.list <- c("_A", "_B", "_C")

####
diss.tt.dat <- read_csv(dissolved.file) %>%
  filter(!str_detect(SampID, "Blk")) %>%
  filter(!str_detect(SampID, "CXBLK")) %>%
  filter(!str_detect(SampID, "Poo")) %>%
  filter(!str_detect(SampID, "MQ")) %>%
  filter(!str_detect(SampID, "Mix")) %>%
  filter(str_detect(SampID, str_c(match.dat$Diss_Samp_ID, collapse="|"))) %>%
  mutate(SampID = str_remove(SampID, "220602_Smp_"),
         SampID = str_remove(SampID, "220623_Smp_"),
         SampID = str_remove(SampID, "231004_Smp_"),
         SampID = str_remove(SampID, "231004_Smp_"),
         SampID = str_remove(SampID, "240102_Smp_"),
         Diss_Samp_ID = str_remove(SampID, "_0nM"),
         Diss_Samp_ID = str_remove(Diss_Samp_ID, str_c(repid.list, collapse="|"))) %>%
  mutate(Diss_Samp_ID = str_replace(Diss_Samp_ID, "_A", "")) %>%
  mutate(Diss_Samp_ID = str_replace(Diss_Samp_ID, "_B", "")) %>%  
  mutate(Diss_Samp_ID = str_replace(Diss_Samp_ID, "_C", "")) %>%
  mutate(Diss_Samp_ID = str_replace(Diss_Samp_ID, "S73", "S7_C3")) %>%
  left_join(., match.dat) %>%
  left_join(., tt.dat) %>%
  filter(Name == Compound) 


region.pal <- c("#00887d", "#014d64", "#01a2d9", "#6794a7", "#76c0c1")



####Flux calculations
flux.dat <- diss.tt.dat %>%
  mutate(comp_flux = (EE.adjust.conc/TT)*24,
         carbon_flux = comp_flux*C) %>%
  filter(!Diss_Samp_ID == "MU12")

TT.sum <- flux.dat %>%
  select(TT, KinExp_ID, Compound, Region) %>%
  unique() %>%
  ungroup() %>%
  mutate(Region = as.factor(Region)) %>%
  mutate(Region = fct_relevel(Region, c("NPTZ", "Gyre", "Equatorial", "Coastal")))

#TT
TT.plot <- ggplot(TT.sum, aes(x = KinExp_ID, y = TT, fill = Region)) +
  facet_wrap(.~Compound, scales = "free_x") +
  theme_test() +
  theme(axis.text.x = element_text(angle = 70, vjust = 0.5)) +
  scale_fill_manual(values = region.pal) +
  geom_col(width = 0.6, linewidth = 0.08, color = "black") +
  scale_y_continuous(expand = c(0,0), limits = c(0,350)) +
  xlab("Experiment") +
  ylab("Turnover Time (hr)")
######

TT.plot
ggsave(TT.plot, dpi = 300, width = 8, height = 4, scale = 1.1, units = "in", file = "Figures/TT_Plot.png")







#Dissolved Conc.
ggplot(flux.dat, aes(x = KinExp_ID, y = EE.adjust.conc)) +
  geom_point() + 
  facet_wrap(.~Compound, scales = "free_x") +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_y_log10()

#Flux
Flux.sum <- flux.dat %>%
  select(carbon_flux, KinExp_ID, Compound, Region) %>%
  group_by(KinExp_ID, Compound, Region) %>%
  summarize(mean.c.flux = mean(carbon_flux),
            sd_carbon_flux = sd(carbon_flux)) 

Flux.sum.2 <- Flux.sum %>%
  ungroup() %>%
  mutate(Region = as.factor(Region)) %>%
  mutate(Region = fct_relevel(Region, c("NPTZ", "Gyre", "Equatorial", "Coastal")))
  #mutate(Region2 = reorder(as.factor(Region), c("NPTZ", "Gyre", "Equatorial", "Coastal")))

library(ggthemes)

flux.fig <- ggplot(Flux.sum.2, aes(x = KinExp_ID, y = mean.c.flux, fill = Region)) +
  geom_point(size = 6, shape = 21) + 
  theme_bw() +
  scale_fill_manual(values = region.pal) +
 # geom_errorbar(aes(ymin = mean.c.flux - sd_carbon_flux, ymax = mean.c.flux + sd_carbon_flux),
 #               width = 0, size = 1) +
  facet_wrap(.~Compound, scales = "free_x") +
  theme(axis.text.x = element_text(angle = 30, vjust = 0.5)) +
  scale_y_log10() +
  xlab("Experiment") +
  ylab("Mean Carbon Flux (nmol/L/day)") 
flux.fig 
ggsave(flux.fig, dpi = 300, width = 8, height = 4, scale = 1.1, units = "in", file = "Figures/Flux_fig.png")

















