

##Culture Figures and Analysis

#install.packages("ggthemes")

library(tidyverse)
library(ggthemes)
library(viridis)

cult.metab.file <- "Intermediates/Culture_Final_Quant_QCed.csv"
cult.meta.data.file <- "Meta_Data/Ingalls_Lab_Data/Culture_Meta_Data.csv"

#####
metab.dat <- read_csv(cult.metab.file) %>%
  rename("Samp_ID" = SampID)
meta.dat <- read_csv(cult.meta.data.file)





### combine metabolite and meta data and calculate per cell and per biovolume concentrations of metabolites
cult.metabs <- left_join(metab.dat, meta.dat) %>%
  mutate(umol_in_vial = uM.in.vial.ave*10^-6*400,
         nmol_per_cell = (umol_in_vial/cells_filt)*1000,
         cell_conc_nM = (umol_in_vial*1000/(cell_volume_on_filter_uL*10^-6)),
         cell_conc_uM = cell_conc_nM/1000,
         cell_conc_mM = cell_conc_uM/1000)

cult.sum <- cult.metabs %>%
  select(Organism, Batch, Compound, cell_conc_mM) %>%
  group_by(Organism, Batch, Compound) %>%
  reframe(mean_conc_mM = mean(cell_conc_mM)) %>%
  filter(!is.na(mean_conc_mM)) %>%
  mutate(Class = as.factor(case_when(Compound %in% c("beta-Alaninebetaine", "Glycine betaine", "Homarine",
                                                   "Trigonelline", "Proline betaine", "Betonicine") ~ "Betaine",
                                   Compound %in% c("Dimethylsulfoniopropionate", "Dimethylsulfonioacetate", "Gonyol") ~ "Sulfonium",
                                   TRUE ~ "Amine_Oxide"))) %>%
  mutate(Class = fct_relevel(Class, c("Betaine", "Sulfonium", "Amine_Oxide")))  %>%
  filter(!Compound == "Trimethylamine N-oxide")


#Visualize culture betaine and sulfonium concentrations:

#stacked bar chart
ggplot(cult.sum, aes(x = Organism, y = mean_conc_mM, fill = Compound)) +
  geom_col() + 
  scale_fill_tableau() +
  facet_wrap(.~Class, nrow = 2, scales = "free_y")

#filled bar chart
ggplot(cult.sum, aes(x = Organism, y = mean_conc_mM, fill = Compound)) +
  geom_col(position = "fill") +
  scale_fill_tableau() +
  facet_wrap(.~Class, nrow = 2)



##Make heatmap:
cult.heat.dat <- cult.metabs %>%
  mutate(Class = as.factor(case_when(Compound %in% c("beta-Alaninebetaine", "Glycine betaine", "Homarine",
                                                     "Trigonelline", "Proline betaine", "Betonicine") ~ "Betaine",
                                     Compound %in% c("Dimethylsulfoniopropionate", "Dimethylsulfonioacetate", "Gonyol") ~ "Sulfonium",
                                     TRUE ~ "Amine_Oxide"))) %>%
  mutate(Class = fct_relevel(Class, c("Betaine", "Sulfonium", "Amine_Oxide"))) %>%
  group_by(Organism, Compound, Batch, Class) %>%
  reframe(ave_vial_conc = mean(uM.in.vial.ave)) %>%
  group_by(Organism, Batch, Class) %>%
  mutate(rel_vial_conc = ave_vial_conc/max(ave_vial_conc))
    
#Concentration in Vial
ggplot(cult.heat.dat, aes(x = Compound, y = Organism, fill = log10(ave_vial_conc))) +
  geom_tile() +
  facet_grid(.~Class, scales = "free_x", space = "free_x") + 
  theme(axis.text.x = element_text(angle = 30, vjust = 0.5, hjust = 0.4)) 

#Relative Concentration
ggplot(cult.heat.dat, aes(x = Compound, y = Organism, fill = log10(rel_vial_conc))) +
  geom_tile() +
  facet_grid(.~Class, scales = "free_x", space = "free_x") + 
  theme(axis.text.x = element_text(angle = 30, vjust = 0.5, hjust = 0.4)) +
  scale_fill_viridis(limits = c(log10(0.0001), log10(1)), oob = scales::squish,
                     option = "G", end = 0.8, begin = 0.1)


betaine.cult.dat <- cult.heat.dat %>%
  filter(Class == "Betaine") %>%
  mutate(Compound = as.factor(Compound),
         Compound = fct_relevel(Compound, c("Proline betaine", "Betonicine", "Trigonelline",
                                            "Homarine", "beta-Alaninebetaine", "Glycine betaine"))) %>%
  group_by(Organism) %>%
  mutate(betaine.rank = rank(1/ave_vial_conc))
         
         
ggplot(betaine.cult.dat, aes(y = Compound, x = rel_vial_conc)) +
  geom_jitter(height = 0.1, width = 0, shape = 21, size = 2.5, aes(fill = Batch)) +
  scale_x_log10() #+
  #facet_grid(.~Class, scales = "free_x", space = "free_x") 

ggplot(betaine.cult.dat, aes(y = Compound, x = betaine.rank)) +
  geom_jitter(height = 0.1, width = 0.1, shape = 21, size = 2.5, aes(fill = Batch)) 

###Make heatmap with hclust by betaine relative conc:
library(vegan)
library(ggdendro)
library(viridis)

metab.clust.dat <- betaine.cult.dat %>%
  select(Compound, Organism, rel_vial_conc) %>%
  unique() %>%
  pivot_wider(id_cols = Organism, names_from = Compound, values_from = rel_vial_conc) %>%
  column_to_rownames(var = "Organism") 

#change NAs to 0s
metab.clust.dat[is.na(metab.clust.dat)] <- 0

#metab.clust.dat.norm <- decostand(metab.clust.dat, method = "pa", margin = 2)

metab.dist <- vegdist(metab.clust.dat, method = "euclidean")

clust.out <- hclust(metab.dist, method = "average")
plot(clust.out)
dend <- as.dendrogram(clust.out)
dend.dat <- dendro_data(dend)
dend.order <- dend.dat$labels %>%
  rename("order" = x, 
         "Organism" = label) %>%
  select(Organism, order)

metab.heatmap.dat <- betaine.cult.dat %>%
  select(Compound, Organism, rel_vial_conc, betaine.rank) %>%
  left_join(dend.order)

ggplot(metab.heatmap.dat, aes(x = Compound, y = reorder(Organism, order), fill = log10(rel_vial_conc))) +
  geom_tile(color = "black", size = 0.05) + 
  theme(axis.text.x = element_text(angle = 30, vjust = 0.5, hjust = 0.4)) +
  scale_fill_viridis(limits = c(log10(0.0001), log10(1)), oob = scales::squish,
                     option = "G", end = 0.8, begin = 0.1)

ggplot(metab.heatmap.dat, aes(x = Compound, y = reorder(Organism, order), fill = betaine.rank)) +
  geom_tile(color = "black", size = 0.05) + 
  theme(axis.text.x = element_text(angle = 30, vjust = 0.5, hjust = 0.4)) +
  scale_fill_viridis(option = "G", end = 0.8, begin = 0.1, direction = -1)














###Make Number of organisms plot
betaine.count.dat <- betaine.cult.dat %>%
  group_by(Batch, Compound) %>%
  reframe(count = n())

ggplot(betaine.count.dat, aes(y = Compound, x = count, fill = Batch)) +
  geom_col(position = "stack", width = 0.5, color = "black", size = 0.05)















































