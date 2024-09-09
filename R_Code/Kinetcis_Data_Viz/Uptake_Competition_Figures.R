

#load packages 
library(tidyverse)
library(ggthemes)
library(ggpubr)
library(ggsci)
library(igraph)
library(broom)
library(rstatix)
library(network)
library(GGally)


source("Functions.R")

### define inputs

#uptake competition rate file:
comp.file <- "Intermediates/Uptake_Competition_Analysis_Output.csv"

#uptake competition literature synthesis 
synth.file <- "Meta_Data/Data_From_Other_Studies/Uptake_Comp_Lit_Synth.csv"





region.pal <- c("#00887d", "#014d64", "#01a2d9", "#6794a7", "#76c0c1")

dat.fig <- read_csv(comp.file)


##Plot Uptake Competition Experiment
comp.fig <- ggplot(dat.fig, aes(x = treatment, y = nM.per.hour.1)) +
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
comp.fig
ggsave(comp.fig, filename = "Figures/Flux_Paper_Figures/Comp_fig.png", 
       scale = 1.2, height = 4, width = 6, units = "in", dpi = 300)



############_____Uptake Network Visualization:

#make inhibitory network

##identify robustly identified inhibitors (n>1)

#read in synth data:
synth.dat <- read_csv(synth.file)

inhib.comps <- synth.dat %>%
  mutate(Inhib_Effect = case_when(Perc_of_Control < 75 & !Signif_Inhib == "No" ~ 1,
                                  TRUE ~ -1)) %>%
  group_by(Compound_inhib) %>%
  summarize(inhib_score = sum(Inhib_Effect)) %>%
  filter(inhib_score > 0)

##make inhbitory dataset 
dat.inhib <- synth.dat %>%
  #  filter(Perc_of_Control < 75 & !Signif_Inhib == "No") %>%
  filter(Compound_inhib %in% inhib.comps$Compound_inhib) %>%
  group_by(Compound_Test, Compound_inhib) %>%
  summarize(count = n()) %>%
  ungroup()

dat.inhib.rev <- dat.inhib %>%
  select(Compound_inhib, Compound_Test)

net.test <- network(dat.inhib.rev, directed = TRUE)

summary(net.test)

ggnet2(net.test, label = TRUE, node.color = "gray",)

x.test <- ggnet2(net.test, 
       arrow.size = 9, arrow.gap = 0.04, label = TRUE,
       layout.exp = 0.3, node.alpha = 0.5, edge.size = 0.3, edge.alpha = 0.6,
       node.color = "#76c0c1")
x.test
ggsave(x.test, filename = "Figures/Flux_Paper_Figures/net_test.png", scale = 1.5,
       height = 4, width = 4, units = "in")

####

#####
test.2 <- ggarrange(NA, comp.fig, NA, NA, NA, NA, NA, x.test, NA,
                    nrow = 3, heights = c(0.45, 0.05, 0.50),
                    ncol = 3, widths = c(0.04, 0.92, 0.04),
                    labels = c(NA, "A", NA, NA, NA, NA, NA, "B", NA))
test.2 
ggsave(test.2, filename = "Figures/Flux_Paper_Figures/comp_net_fig.png", scale = 1.2,
height = 7, width = 6, units = "in", bg = "white")






####Analyze Uptake Rate Data and Create barchart 

##Subset data 
#dat.comp <- read_csv(comp.file) %>%
#  filter(exp == "UCH1") %>%
#  filter(!str_detect(SampID, "Blk")) %>%
##  select(SampID, cruise, exp, treatment, rep, Fragment_mz, nM.per.hour.1) %>%
#  filter(Fragment_mz == 78.0399) %>%
#  filter(!str_detect(SampID, "TN397_UCH1_Hom")) %>%
#  group_by(cruise, treatment) %>%
#  mutate(mean.rate = mean(nM.per.hour.1)) %>%
#  group_by(cruise, exp) %>%
#  mutate(perc.of.control = mean.rate/mean.rate[treatment == "Gluc"])



####Run test in separate groups, arrange groups by "control" type compound, 
# run separate fdr correction on only comparisons we are interested in 
# join together
# dat.t <- dat.comp %>%
#   arrange(desc(nM.per.hour.1)) %>%
#   group_by(cruise) %>%
#   t_test(nM.per.hour.1~treatment, ref.group = "Gluc")
# 
# 
# dat.gluc <- data.frame(
#   cruise = c("RC104", "TN397"),
#   treatment = c("Gluc", "Gluc"),
#   signif = c("not significant", "not significant"))
# 
# dat.t.res <- dat.t %>%
#   select(cruise, group2, p.adj) %>%
#   mutate(signif = case_when(p.adj > 0.1 ~ "not significant",
#                             p.adj <= 0.1 ~ "significant")) %>%
#   rename("treatment" = group2) %>%
#   select(-p.adj) %>%
#   rbind(dat.gluc)
# 
# #combine statistical results with figure data:
# dat.fig <- dat.comp %>%
#   left_join(., dat.t.res)



# 
# ###define nodes
# nodes.b <- dat.inhib %>%
#   select(Compound_Test) %>%
#   unique() %>%
#   rename("node" = Compound_Test)
# 
# nodes.a <- dat.inhib %>%
#   select(Compound_inhib) %>%
#   unique() %>%
#   rename("node" = Compound_inhib)
# 
# nodes.all <- full_join(nodes.a, nodes.b)
# 
# ###define edges
# edges <- dat.inhib %>%
#   select(Compound_inhib, Compound_Test, count)
# 
# net <- graph_from_data_frame(d=edges, vertices = nodes.all, directed = T)
# 
# #Make edge withs proportional to number of experiments
# E(net)$width <- (E(net)$count^1)
# 
# ###plot network
# set.seed = 124
# 
# net.plot <- plot(net, vertex.size = 10, vertex.label.dist = 4, edge.arrow.size = 0.8,
#      vertex.color = "#00887d", edge.color = "gray65", vertex.label.family = "Helvetica",
#      vertex.label.color = "black")
# net.plot
# 
# install.packages("GGally")
#library(GGally)


















































































