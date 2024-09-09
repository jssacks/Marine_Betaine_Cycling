
#This script produces main text figure with kinetics curves and flux measurments

#load packages 
library(tidyverse)
library(ggthemes)
library(ggpubr)
library(ggsci)
library(ggbreak)
library(png)
library(cowplot)
library(ggrepel)

###Define inputs:
####Define inputs:
flux.file <- "Intermediates/Flux_Analysis_Output.csv"
kin.file <- "Intermediates/Kin_Analysis_Output.csv"


#load in kinetics data:
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
  mutate(Region = fct_relevel(Region, c("NPTZ", "NPSG", "Equatorial", "Puget Sound"))) %>%
  mutate(grid.val = "b") 

grid.add.1 <- kin.curve.dat %>%
  filter(nM_per_hour_noS < 0.8)

grid.add.2 <- kin.curve.dat %>%
  filter(Cruise_Exp == "RC078_UKG1") %>%
  mutate(grid.val = "a") %>%
  filter(nM_per_hour_noS > 1)

kin.fig.dat <- rbind(grid.add.1, grid.add.2)

###Plot kinetics Curves:
region.pal <- c("#00887d", "#014d64", "#01a2d9", "#6794a7")

kin.curves.plot <- ggplot(kin.fig.dat, aes(x = treatment_nM, y = nM_per_hour_noS, color = Region, group = Cruise_Exp)) +
  geom_line(size = 1, alpha = 0.8) +
  #facet_wrap(Compound~.) +
  scale_color_manual(values = region.pal) +
  theme_test() +
  facet_grid(grid.val~Compound, scales = "free", space = "free") +
  #  scale_y_sqrt() +
  
  
  # scale_y_cut(breaks = c(1, 3)) +
  # scale_y
  # scale_y_break(c(0.85, 1.0), scales = "0.40", 
  #                space = 0.06, ticklabels = c(1.0, 2.0, 3.0, 4.0),
  #               expand = c(0, NA, NA, 1000)) +
  ylab(expression(Uptake~Rate~(nmol~L^-1~day^-1))) +
  xlab("Spike Concentration (nM)") +
  theme(strip.background.y = element_blank(),
        strip.text.y = element_blank())

kin.curves.plot



########Try inset plot:
k.d.a <- kin.curve.dat %>%
  filter(treatment_nM <= 1000) %>%
  filter(nM_per_hour_noS < 0.9)

k.d.b <- kin.curve.dat %>%
  filter(Cruise_Exp == "RC078_UKG1")


label.dat <- k.d.a %>%
  filter(treatment_nM == 1000) %>%
  mutate(x = 1000,
         y = nM_per_hour_noS) %>%
  select(Cruise_Exp, Region, x, y, Compound) 





kin.plot.A <- ggplot(k.d.a, aes(x = treatment_nM, y = nM_per_hour_noS, color = Region, group = Cruise_Exp)) +
  geom_line(size = 1, alpha = 0.8) +
  #facet_wrap(Compound~.) +
  scale_color_manual(values = region.pal) +
  theme_test() +
  facet_grid(.~Compound) +
  ylab(expression(Uptake~Rate~(nmol~L^-1~day^-1))) +
  xlab("Spike Concentration (nM)") +
  theme(#strip.background.x = element_blank(),
    strip.text.x = element_text(face = "bold")
  ) +
  xlim(0, 1400) +
  scale_y_continuous(expand = c(0, NA, NA, NA), limits = c(0, 0.9)) +
  geom_text_repel(data = label.dat, aes(x = x, y = y, label = Cruise_Exp),
                  size = 2.5,
                  xlim = c(1000, 1500),
                  min.segment.length = 0.01,
                  color = "black",
                  segment.color = "black",
                  segment.size = 0.25)

#geom_label(data = label.dat, aes(x = x, y = y, label = Cruise_Exp), size = 1)
kin.plot.A

kin.plot.B  <- ggplot(k.d.b, aes(x = treatment_nM, y = nM_per_hour_noS)) +
  geom_line(size = 1, alpha = 0.8, color = "#6794a7") +
  theme_test() +
  theme(axis.title = element_blank()) +
  #  facet_grid(.~Compound) +
  ylab(expression(Uptake~Rate~(nmol~L^-1~day^-1))) +
  xlab("Spike Concentration (nM)") +
  theme(panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA)) +
  annotate("text", x = 4000, y = 1.5, label = "RC078_UKG1", size = 2.5)

kin.plot.B

kin.plot.comb.2 <- ggdraw(kin.plot.A) +
  draw_plot(kin.plot.B,
            x = .2, y = .7, .2, .2, scale = 1.2) 
kin.plot.comb.2

#save plot
ggsave(kin.plot.comb.2, file = "Figures/Flux_Paper_Figures/Kinetics_Curves_inset.pdf",
       height = 3.5, width = 6, dpi = 300, units = "in", scale = 1.3, bg = "white")
