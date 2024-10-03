
#This script produces a combined uptake kinetics experiment and uptake competition experiment figure


#load packages 
library(tidyverse)
library(ggthemes)
library(ggpubr)
#library(ggsci)
#library(ggbreak)
#library(png)
library(cowplot)
library(ggrepel)

###Define inputs:

#kinetics data
kin.file <- "Intermediates/Kin_Analysis_Output.csv"

#uptake competition rate file:
comp.file <- "Intermediates/Uptake_Competition_Analysis_Output.csv"



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
                  segment.size = 0.25) +
  theme(legend.position = "bottom")

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
            x = .19, y = .75, .15, .15, scale = 1.4) 
kin.plot.comb.2



####Make uptake competition plot:
#region.pal <- c("#00887d", "#014d64", "#01a2d9", "#6794a7", "#76c0c1")

dat.fig <- read_csv(comp.file)


##Plot Uptake Competition Experiment
comp.plot <- ggplot(dat.fig, aes(x = treatment, y = nM.per.hour.1)) +
  geom_col(aes(x = reorder(treatment, order), y = mean.rate, fill = signif), 
           color = "black", size = 0.1, width = 0.6, position = "dodge") +
  scale_fill_manual(values = c("lightgray", "#76c0c1",  "#00887d")) +
  geom_jitter(width = 0.09, shape = 21, fill = "white", size = 2) +
  theme_test() +
  facet_grid(.~cruise, scales = "free_x", space = "free") +
  xlab("Treatment") +
  ylab("Uptake Rate (nM/hr)") +
  labs(fill = "") + 
  theme(legend.position = "bottom") + 
  scale_y_continuous(expand = c(0, NA, NA, 1000), limits = c(0,0.11)) 
comp.plot


###Combine plots together:
exp.fig <- ggarrange(NA, NA, NA,
                     NA, comp.plot, NA,
                     NA, NA, NA,
                     NA, kin.plot.comb.2, NA,
                     NA, NA, NA,
                     nrow = 5, ncol = 3,
                     heights = c(0.02, 0.38, 0.04, 0.54, 0.02),
                     widths = c(0.04, 0.92, 0.04),
                     labels = c(NA, NA, NA,
                                NA, "A", NA,
                                NA, NA, NA,
                                NA, "B", NA,
                                NA, NA, NA))
exp.fig

ggsave(exp.fig, filename = "Figures/Flux_Paper_Figures/Kinetics_Competition_Experiments_Figure.png",
       dpi = 600, units = "in", scale = 1.3, bg = "white", height = 7, width = 6)



































