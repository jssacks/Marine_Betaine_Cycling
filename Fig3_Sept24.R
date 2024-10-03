

#packages:
library(tidyverse)
library(ggpubr)
library(cowplot)



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
  left_join(., region.dat) %>%
  mutate(norm_kt_diff = (mean_ks - Mean.Diss.Conc.nM)/(mean_ks)*100,
         perc_diss_pool = (Mean.Diss.Conc.nM/metab.tot.nM)*100) %>%
  mutate(Compound = as.factor(Compound),
         Region = as.factor(Region))





# Ks vs. Compound Dissolved Concentration of Homarine
all.dat.hom <- all.dat %>%
  filter(Compound == "Homarine")

plot.a <- ggplot(all.dat.hom, aes(x = Mean.Diss.Conc.nM, y = mean_ks)) +
  geom_abline(slope = 1, linetype = "dashed") +
  geom_smooth(alpha = 0.2, method = "lm", color = "gray20") +
  geom_errorbar(aes(ymin = mean_ks-sd_ks, ymax = mean_ks+sd_ks), width = 0.03) +
  geom_errorbarh(aes(xmin = Mean.Diss.Conc.nM-SD.Diss.Conc.nM, xmax = Mean.Diss.Conc.nM+SD.Diss.Conc.nM), height = 0.03) +
  geom_point(aes(shape = factor(Compound, levels = c("GBT", "Homarine")), fill = factor(Region, levels = c("NPTZ", "NPSG", "Equatorial", "Puget Sound"))), size = 3.5) +
  scale_shape_manual(values = c("GBT" = 21, "Homarine" = 24), drop = FALSE) +
  scale_fill_manual(values = c("NPTZ" =  "#00887d",
                               "NPSG" =  "#014d64",
                               "Equatorial" = "#01a2d9",
                               "Puget Sound" = "#6794a7"), drop = FALSE) +
  scale_x_log10() +
  scale_y_log10(limits = c(1, 1000)) +
  theme_test()  +
  guides(fill = guide_legend(override.aes = list(shape = 21)),
         shape = guide_legend(override.aes = list(fill = "black"))) +
  xlab("Dissolved Concentration (nM)") +
  ylab(expression(K[t]~(nM))) +
  annotate("text", x = 0.4, y = 5, 
           label = "p = 0.09",
           size = 3)  
  
plot.a


# Ks vs. Compound Dissolved Concentration of Homarine
all.dat.gbt <- all.dat %>%
  filter(Compound == "GBT")

plot.b <- ggplot(all.dat.gbt, aes(x = Mean.Diss.Conc.nM, y = mean_ks)) +
  geom_abline(slope = 1, linetype = "dashed") +
  geom_smooth(alpha = 0.2, method = "lm", color = "gray20") +
  geom_errorbar(aes(ymin = mean_ks-sd_ks, ymax = mean_ks+sd_ks), width = 0.03) +
  geom_errorbarh(aes(xmin = Mean.Diss.Conc.nM-SD.Diss.Conc.nM, xmax = Mean.Diss.Conc.nM+SD.Diss.Conc.nM), height = 0.03) +
  geom_point(aes(shape = Compound, fill = Region), size = 3.5) +
  scale_shape_manual(values = c("GBT" = 21, "Homarine" = 24), drop = FALSE) +
  scale_fill_manual(values = c("NPTZ" =  "#00887d",
                               "NPSG" =  "#014d64",
                               "Equatorial" = "#01a2d9",
                               "Puget Sound" = "#6794a7"), drop = FALSE) +
  scale_x_log10() +
  scale_y_log10(limits = c(1, 1000)) +
  theme_test()  +
  guides(fill = guide_legend(override.aes = list(shape = 21)),
         shape = guide_legend(override.aes = list(fill = "black"))) +
  xlab("Dissolved Concentration (nM)") +
  ylab(expression(K[t]~(nM))) +
  annotate("text", x = 1.1, y = 400, 
           label = "p = 0.01",
           size = 3)  

plot.b



######Relative Kt difference plotted against % of total betaine + sulfonium pool
plot.c <- ggplot(all.dat, aes(x = perc_diss_pool, y = norm_kt_diff)) + 
  geom_smooth(alpha = 0.2, method = "lm", color = "gray20") +
  geom_point(aes(shape = Compound, fill = Region), size = 3.5) +
  scale_shape_manual(values = c(21, 24), drop = FALSE) +
  scale_fill_manual(values = region.pal, drop = FALSE)  +
  theme_test()  +
  guides(fill = guide_legend(override.aes = list(shape = 21)),
         shape = guide_legend(override.aes = list(fill = "black"))) +
  xlab("Percentage of Dissolved Pool (%)") +
  ylab("Normalized Difference (% of Kt)") +
  annotate("text", x = 10, y = 80, 
           label = "p = 0.0007",
           size = 3) 

plot.c



# Flux vs. Particulate Concentration
plot.d <- ggplot(all.dat, aes(x = mean.part.conc.nM, y = mm_flux_nM_day)) +
  geom_abline(slope = 1, linetype = "dashed") +
  geom_smooth(alpha = 0.2, method = "lm", color = "gray20") +
  geom_errorbar(aes(ymin = mm_flux_nM_day-mm_flux_sd_nM_day, ymax = mm_flux_nM_day+mm_flux_sd_nM_day), width = 0.04) +
  geom_errorbarh(aes(xmin = mean.part.conc.nM-sd.part.conc.nM, xmax = mean.part.conc.nM+sd.part.conc.nM), height = 0.04) +
  geom_point(aes(shape = Compound, fill = Region), size = 3.5) +
  scale_shape_manual(values = c(21, 24)) +
  scale_fill_manual(values = region.pal) +
  scale_x_log10() +
  scale_y_log10() +
  theme_test()  +
  guides(fill = guide_legend(override.aes = list(shape = 21)),
         shape = guide_legend(override.aes = list(fill = "black"))) +
  xlab("Particulage GBT or Homarine Concentration (nM)") +
  ylab(expression(Flux~(nmol~L^-1~day^-1))) +
  annotate("text", x = 0.05, y = 20, 
           label = "p = 0.0004",
           size = 3) +
  theme(legend.position = "top") 
plot.d

g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

all.legend <- g_legend(plot.d)


plot.a.nl <- plot.a +
  theme(legend.position = "none")
plot.b.nl <- plot.b +
  theme(legend.position = "none")
plot.c.nl <- plot.c +
  theme(legend.position = "none")
plot.d.nl <- plot.d +
  theme(legend.position = "none")

###combine all four plots:
plot.all <- ggarrange(ggarrange(NA, NA, NA, NA, NA,
                      NA, plot.a.nl, NA, plot.b.nl, NA,
                      NA, NA, NA, NA, NA,
                      NA, plot.c.nl, NA, plot.d.nl, NA,
                      NA, NA, NA, NA, NA,
                     # common.legend = TRUE,
                    #  legend = "bottom",
                      nrow = 5, ncol = 5,
                      widths = c(0.02, 0.43, 0.02, 0.43, 0.02),
                      heights = c(0.02, 0.43, 0.02, 0.43, 0.02),
                      labels = c(NA, NA, NA, NA, NA,
                                 NA, "A", NA, "B", NA,
                                 NA, NA, NA, NA, NA,
                                 NA, "C", NA, "D", NA,
                                 NA, NA, NA, NA, NA)),
                    all.legend, nrow = 2, ncol = 1,
                    heights = c(0.9, 0.1))
plot.all

ggsave(plot.all, filename = "Figures/Flux_Paper_Figures/Figure3.png",
       scale = 1.25, units = "in", dpi = 600,
       width = 7, height = 6.5, bg = "white")




# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# #get legend
# legend <- get_legend(plot.a)
# 
# #remove legend from all plots
# plot.a2 <- plot.a +
#   theme(legend.position = "none")
# 
# plot.b2 <- plot.b +
#   theme(legend.position = "none")
# 
# plot.c2 <- plot.c +
#   theme(legend.position = "none")
# 
# plot.all <- ggarrange(NA, NA, NA, NA, NA,
#                       NA, plot.a2, NA, plot.b2, NA,
#                       NA, NA, NA, NA, NA,
#                       NA, plot.c2, NA, legend, NA,
#                       NA, NA, NA, NA, NA,
#                       nrow = 5, ncol = 5,
#                       widths = c(0.02, 0.43, 0.02, 0.43, 0.02),
#                       heights = c(0.02, 0.43, 0.02, 0.43, 0.02),
#                       labels = c(NA, NA, NA, NA, NA,
#                                  NA, "A", NA, "B", NA,
#                                  NA, NA, NA, NA, NA,
#                                  NA, "C", NA, NA, NA,
#                                  NA, NA, NA, NA, NA))
# 
# plot.all
# ggsave(plot.all, filename = "Figures/Flux_Paper_Figures/KT_DissConc_Flux_PartConc_LM_Figure.png",
#        scale = 1.2, units = "in", dpi = 600,
#        width = 7, height = 6, bg = "white")