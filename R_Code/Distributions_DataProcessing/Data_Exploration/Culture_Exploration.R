


library(tidyverse)
library(viridis)

cult.dat <- read_csv("Intermediates/Culture_Quant_Output.csv")
dat.info <- read_csv("Meta_Data/Ingalls_Lab_Data/Culture_Meta_Data.csv")
dat.blk.key <- read_csv("Meta_Data/Ingalls_Lab_Data/Culture_Blk_Match_Key.csv")

dat.sml <- cult.dat %>%
  rename("Samp_ID" = SampID) %>%
  left_join(., dat.info) %>%
  filter(!is.na(Organism)) %>%
  group_by(Organism, Name, Batch) %>%
  summarize(mean.conc = mean(umol.in.vial.ave)) %>%
  ungroup()

dat.blks <- cult.dat %>%
  filter(str_detect(SampID, "Blk")) %>%
  rename("Blk_ID" = SampID) %>%
  left_join(., dat.blk.key) %>%
  group_by(Name, Organism) %>%
  mutate(mean.blk = mean(umol.in.vial.ave),
         blk.threshold = 3*mean.blk) %>%
  select(Name, Organism, blk.threshold) %>%
  unique()

dat.blk.filt <- dat.sml %>%
  ungroup() %>%
  left_join(., dat.blks) %>%
  mutate(blk.flag = case_when(mean.conc <= blk.threshold ~ "Flag",
                              TRUE ~ NA)) %>%
  filter(is.na(blk.flag))

ggplot(dat.sml, aes(x = reorder(Organism, Batch), y = mean.conc, fill = Batch)) +
  geom_col() +
  facet_wrap(.~Name, scales = "free")

dat.org.rel <- dat.blk.filt %>%
  filter(!Name %in% c("Gonyol", 
                      "Dimethylsulfoniopropionate", 
                      "Dimethylsulfonioacetate",
                      "Trimethylamine N-oxide")) %>%
  group_by(Organism) %>%
  mutate(Rel.Conc = mean.conc/max(mean.conc))

ggplot(dat.org.rel, aes(x = Name, y = reorder(Organism, Batch), fill = log10(Rel.Conc))) +
  geom_tile() +
  scale_fill_viridis()

ggplot(dat.org.rel, aes(x = Name, y = mean.conc)) +
  geom_boxplot() +
  geom_jitter(aes(color = Batch), width = 0.1) +
  scale_y_log10()

ggplot(dat.org.rel, aes(x = Name, y = mean.conc)) +
  geom_boxplot() +
  facet_wrap(.~Batch) +
  geom_jitter(aes(color = Batch), width = 0.1) +
  scale_y_log10()


####QCed data
qced.dat <- read_csv("Intermediates/Culture_betaine_QCdat_samplesremoved.csv") %>%
  group_by(Compound, Organism, Batch) %>%
  summarize(mean.area = mean(Area)) %>%
  ungroup()





####Do hclust:
#perform clustering:
library(vegan)
library(ggdendro)
library(viridis)

metab.clust.dat <- qced.dat %>%
  select(Compound, Organism, mean.area) %>%
  unique() %>%
  pivot_wider(id_cols = Organism, names_from = Compound, values_from = mean.area) %>%
  column_to_rownames(var = "Organism") 

#change NAs to 0s
metab.clust.dat[is.na(metab.clust.dat)] <- 0

metab.clust.dat.norm <- decostand(metab.clust.dat, method = "pa", margin = 2)

metab.dist <- vegdist(metab.clust.dat.norm, method = "euclidean")

clust.out <- hclust(metab.dist, method = "average")
plot(clust.out)
dend <- as.dendrogram(clust.out)
dend.dat <- dendro_data(dend)
dend.order <- dend.dat$labels %>%
  rename("order" = x, 
         "Organism" = label) %>%
  select(Organism, order)

metab.heatmap.dat <- qced.dat %>%
  select(Compound, Organism, mean.area) %>%
  left_join(dend.order)



ggplot(metab.heatmap.dat, aes(x = Compound, y = reorder(Organism, order), fill = log10(mean.area))) +
  geom_tile(color = "black") +
  scale_fill_viridis() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
  theme_test() 



