#This script produces main text figure with kinetics curves and flux measurments

#load packages 
library(tidyverse)
library(ggthemes)
library(ggpubr)

###Define inputs:
####Define inputs:
kin.file <- "Intermediates/Kin_Analysis_Output.csv"
exp.file <- "Intermediates/All_Kin_Exp_Dat.csv"

region.pal <- c("#00887d", "#014d64", "#01a2d9", "#6794a7")

#load in kinetics data:
exp.dat <- read_csv(exp.file)
kin.dat <- read_csv(kin.file)

##Generate Experimental Kinetics Curves based on fit MM equations
##Create a dataframe to plot the nonlinear equations with

spike.vals <- data.frame(treatment_nM = seq(1, 6000, by = 1)) %>%
  mutate(join.val = "x")

kin.curve.dat <- kin.dat %>%
  select(Cruise, exp, mean_ks, mean_vmax, mean_s) %>% 
  unite(c("Cruise", "exp"), col = "Cruise_Exp", remove = FALSE) %>%
  mutate(join.val = "x") %>%
  left_join(., spike.vals, relationship = "many-to-many") %>%
  mutate(nM_per_hour = mean_vmax*(treatment_nM+mean_s)/(mean_ks+treatment_nM+mean_s)) %>%
  mutate(nM_per_hour_noS = (mean_vmax*(treatment_nM))/(mean_ks+treatment_nM)) %>%
  mutate(Compound = case_when(str_detect(Cruise_Exp, "UKG") ~ "Glycine betaine",
                              str_detect(Cruise_Exp, "UKH") ~ "Homarine",
                              str_detect(Cruise_Exp, "UCH") ~ "Homarine")) %>%
  mutate(Region = case_when(str_detect(Cruise_Exp, "RC104") ~ "Puget Sound",
                            str_detect(Cruise_Exp, "RC078") ~ "Puget Sound",
                            str_detect(Cruise_Exp, "KM1906") ~ "NPTZ",
                            str_detect(Cruise_Exp, "TN412") ~ "NPSG",
                            str_detect(Cruise_Exp, "TN397_UKH1") ~ "NPSG",
                            str_detect(Cruise_Exp, "TN397_UKH1") ~ "NPSG",
                            str_detect(Cruise_Exp, "TN397_UKH2") ~ "NPSG",
                            str_detect(Cruise_Exp, "TN397_UKH3") ~ "Equatorial",
                            str_detect(Cruise_Exp, "TN397_UKH4") ~ "Equatorial",
                            str_detect(Cruise_Exp, "TN397_UKH5") ~ "Equatorial")) %>%
  mutate(Region = as.factor(Region)) %>%
  mutate(Region = fct_relevel(Region, c("NPTZ", "NPSG", "Equatorial", "Puget Sound")))








#####Pull in and Organize Data from experiments:
exp.dat.tidy <- exp.dat %>%
  unite(c("cruise", "exp"), col = "Cruise_Exp", remove = FALSE) %>%
  filter(Fragment_mz %in% c(97.0954, 62.8977))


exp.max.spike <- exp.dat.tidy%>%
  group_by(Cruise_Exp) %>%
  mutate(max_spike = max(treatment_conc)) %>%
  select(Cruise_Exp, max_spike) %>%
  unique()


###Organize MM line data:
kin.curve.dat.to.plot <- kin.curve.dat %>%
  left_join(., exp.max.spike) %>%
  filter(treatment_nM <= max_spike)




all.exp.plot <- ggplot(exp.dat.tidy, aes(x = treatment_conc, y = nM.per.hour.1)) +
  facet_wrap(.~Cruise_Exp, scales = "free") +
  theme_bw() +
  geom_line(data = kin.curve.dat.to.plot,
            aes(x = treatment_nM, y = nM_per_hour_noS, color = Region),
            size = 1) +
  geom_point(shape = 21, fill = "gray 90", size = 2) +
  scale_color_manual(values = region.pal) +
  theme(legend.position = "bottom") +
  xlab("Spike Concentration (nM)") +
  ylab(expression(Uptake~Rate~(nmol~L^-1~day^-1)))
all.exp.plot


ggsave(all.exp.plot, filename = "Figures/Flux_Paper_Figures/All_Kin_Exp_Supplemental_Figure.png",
       height = 5, width = 6, dpi = 300, units = "in", scale = 1.3, bg = "white")
