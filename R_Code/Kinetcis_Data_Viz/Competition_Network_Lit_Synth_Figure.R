

#Uptake Competition Network File:


#load packages 
library(tidyverse)
library(ggthemes)
library(ggsci)
library(igraph)
library(broom)
library(rstatix)
library(network)
library(GGally)


source("Functions.R")

### define inputs

#uptake competition literature synthesis 
synth.file <- "Meta_Data/Data_From_Other_Studies/Uptake_Comp_Lit_Synth.csv"




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
ggsave(x.test, filename = "Figures/Flux_Paper_Figures/Uptake_Competition_Network.png", scale = 1.3,
       height = 4, width = 5, units = "in")
