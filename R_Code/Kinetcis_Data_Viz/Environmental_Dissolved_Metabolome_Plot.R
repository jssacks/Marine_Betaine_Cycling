


#load packages 
library(tidyverse)
library(ggthemes)
library(ggpubr)
library(ggsci)




#define inputs 
env.file <- "Intermediates/Kinetics_Exp_Enviro_Metabolome_Dat.csv"


#Load data
env.dat <- read_csv(env.file) %>%
  unite(c("Cruise", "exp"), col = "Cruise_Exp")

#Reroder compounds and samples:
ord.comp <- env.dat %>%
  group_by(Compound) %>%
  reframe(mean.conc = mean(Mean.Diss.Conc.nM)) %>%
  mutate(order.comp = rank(-mean.conc)) %>%
  select(Compound, order.comp)

ord.samp <- env.dat %>%
  mutate(exp.comp = case_when(str_detect(Cruise_Exp, "UKG") ~ "Glycine betaine",
                              str_detect(Cruise_Exp, "UKH") ~ "Homarine",
                              str_detect(Cruise_Exp, "UCH") ~ "Homarine")) %>%
  select(Cruise_Exp, exp.comp) %>%
  unique()

env.dat.ord <- env.dat %>%
  left_join(., ord.comp) %>%
  left_join(., ord.samp)

####
plot.a <- ggplot(env.dat.ord, aes(x = reorder(Cruise_Exp, -Mean.Diss.Conc.nM), y = Mean.Diss.Conc.nM, fill = reorder(Compound, order.comp))) +
  geom_col(color = "black", size = 0.05, width = 0.7) +
  scale_fill_npg() +
  facet_grid(.~exp.comp, scales = "free_x", space = "free_x") +
  theme_test() +
  xlab("Experiment") +
  ylab("Dissolved Concentration (nM)") +
  labs(fill = "Compound") +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.8, hjust = 0.8)) +
  scale_y_continuous(expand = c(0, NA, NA, 1000), limits = c(0,160)) 
plot.a

plot.b <- ggplot(env.dat.ord, aes(x = reorder(Cruise_Exp, -Mean.Diss.Conc.nM), y = Mean.Diss.Conc.nM, fill = reorder(Compound, order.comp))) +
  geom_col(color = "black", size = 0.05, position = "fill", width = 0.6) +
  scale_fill_npg() +
  facet_grid(.~exp.comp, scales = "free_x", space = "free_x") +
  theme_test() +
  xlab("Experiment") +
  ylab("Relative Concentration (Mole %)") +
  labs(fill = "Compound") +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.8, hjust = 0.8)) 
plot.b

plot.ab <- ggarrange(plot.a, NA, plot.b,
                     nrow= 3, ncol = 2,
                     heights = c(0.475, 0.475, 0.05, 0.05, 0.475, 0.475),
                     widths = c(0.95, 0.05),
                     common.legend = TRUE,
                     legend = "right",
                     align = "h")

plot.ab

ggsave(plot.ab, filename = "Figures/Flux_Paper_Figures/Environmental_Dissolved_Metabolomes.png",
       height = 8, width = 8, dpi = 300, units = "in", scale = 1.15, bg = "white")

















































































