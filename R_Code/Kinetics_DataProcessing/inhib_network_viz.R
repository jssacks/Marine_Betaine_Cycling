

library(tidyverse)
library(igraph)


#data
synth.dat <- read_csv("Meta_Data/Data_From_Other_Studies/Uptake_Comp_Lit_Synth.csv")



#make inhibitory network

##identify robustly identified inhibitors (n>1)
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



###define nodes
nodes.b <- dat.inhib %>%
  select(Compound_Test) %>%
  unique() %>%
  rename("node" = Compound_Test)

nodes.a <- dat.inhib %>%
  select(Compound_inhib) %>%
  unique() %>%
  rename("node" = Compound_inhib)

nodes.all <- full_join(nodes.a, nodes.b)

###define edges
edges <- dat.inhib %>%
  select(Compound_inhib, Compound_Test, count)

net <- graph_from_data_frame(d=edges, vertices = nodes.all, directed = T)

#Make edge withs proportional to number of experiments
E(net)$width <- (E(net)$count^1)

###plot network
plot(net, vertex.size = 10, vertex.label.dist = 4, edge.arrow.size = 0.8)




########Make noninhibitory network:
#make inhibitory network

##identify robustly identified inhibitors (n>1)
noninhib.comps <- synth.dat %>%
  filter(Inhib == "No") %>%
  group_by(Compound_inhib) %>%
  summarize(count = n()) %>%
  filter(count > 1)

##make inhbitory dataset 
dat.noninhib <- synth.dat %>%
  filter(Inhib == "No") %>%
  filter(Compound_inhib %in% noninhib.comps$Compound_inhib) %>%
  group_by(Compound_Test, Compound_inhib) %>%
  summarize(count = n()) %>%
  ungroup()



###define nodes
n.nodes.a <- dat.noninhib %>%
  select(Compound_Test) %>%
  unique() %>%
  rename("node" = Compound_Test)

n.nodes.b <- dat.noninhib %>%
  select(Compound_inhib) %>%
  unique() %>%
  rename("node" = Compound_inhib)

n.nodes.all <- full_join(n.nodes.a, n.nodes.b)

###define edges
n.edges <- dat.noninhib

n.net <- graph_from_data_frame(d=n.edges, vertices = n.nodes.all, directed = F)

#Make edge withs proportional to number of experiments
E(n.net)$width <- (E(n.net)$count^1)*2

###plot network
plot(n.net, vertex.size = 8, vertex.label.dist = 4, edge.color = "darkred")
































