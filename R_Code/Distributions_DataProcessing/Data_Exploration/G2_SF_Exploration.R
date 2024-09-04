
library(tidyverse)
library(viridis)


#g2.dat <- read_csv("Intermediates/G2_SizeFractionated_HILIC_Pos_BMISed_dat.csv")
g2.metadat <- read_csv("Meta_Data/Ingalls_Lab_Data/G2_SF_Sample_MetaData.csv")
g2.final.dat <- read_csv("Intermediates/G2_Final_Areas_QCed.csv")




g2.blks <- read_csv("Intermediates/g2_betaine_data_raw.csv") %>%
  filter(str_detect(Rep, "Blk")) %>%
  mutate(Blk_type = case_when(str_detect(Rep, "MQ") ~ "MQ",
                              TRUE ~ "FSW")) %>%
  group_by(Compound, Blk_type) %>%
  summarize(ave.blk = mean(Area)) %>%
  rename("MF" = Compound)






#Summarizing, normalizing to vol filtered, and calculating percentage of signal in small size fraction
g2.all.dat <- left_join(g2.dat, g2.metadat) %>%
  filter(!is.na(Station)) %>%
  mutate(Vol.norm.area = Adjusted_Area/Vol_L) %>%
  group_by(MF, Station, Lat, Size_Fraction) %>%
  summarize(Mean.Area = mean(Vol.norm.area)) %>%
  group_by(MF, Station) %>%
  mutate(Sml_Frac_Perc = Mean.Area[Size_Fraction == "S"]/Mean.Area[Size_Fraction == "F"])

###Combining ___ and ___ data to determine if compounds were detected or not:
g2.blk.qc <- full_join(g2.dat, g2.blks) %>%
  left_join(., g2.metadat) %>%
  ungroup() %>%
  filter(!is.na(Station)) %>%
  mutate(blk.ratio = Adjusted_Area/ave.blk)




####
g2.f.meta.dat <- g2.final.dat %>%
 # rename("SampID" = Rep) %>%
  left_join(., g2.metadat %>%
              rename("Rep" = SampID)) %>%
  mutate(Vol.norm.area = Adjusted_Area/Vol_L) %>%
  group_by(Compound, Station, Lat, Size_Fraction) %>%
  reframe(Mean.Area = mean(Vol.norm.area)) %>%
  group_by(Compound, Station, Lat) %>%
  mutate(Mean.Area = case_when(Size_Fraction == "F" ~ Mean.Area-Mean.Area[Size_Fraction == "S"],
                               TRUE ~ Mean.Area)) %>%
  mutate(Size_Fraction = str_replace(Size_Fraction, "F", "L")) %>%
  mutate(Small_Frac_Perc = Mean.Area[Size_Fraction == "S"]/(Mean.Area[Size_Fraction == "S"] + Mean.Area[Size_Fraction == "L"])) %>%
  ungroup()


##Make stacked bar chart:
ggplot(g2.f.meta.dat, aes(x = as.factor(Lat), y = Mean.Area, fill = Size_Fraction)) +
  geom_col(width = 0.5) +
  facet_wrap(.~Compound, scales = "free") +
  xlab("Latitude") +
  ylab("Peak Area")


###Perform hclust to organize samples for heatmap:

#prep dataset
metab.clust.dat <- g2.f.meta.dat %>%
  select(Compound, Lat, Small_Frac_Perc) %>%
  unique() %>%
  pivot_wider(id_cols = Compound, names_from = Lat, values_from = Small_Frac_Perc) %>%
  column_to_rownames(var = "Compound")

#perform clustering:
library(vegan)
library(ggdendro)
library(viridis)

metab.dist <- vegdist(metab.clust.dat, method = "euclidean")

clust.out <- hclust(metab.dist, method = "average")
dend <- as.dendrogram(clust.out)
dend.dat <- dendro_data(dend)
dend.order <- dend.dat$labels %>%
  rename("order" = x, 
         "Compound" = label) %>%
  select(Compound, order)

metab.heatmap.dat <- g2.f.meta.dat %>%
  select(Compound, Lat, Small_Frac_Perc) %>%
  left_join(dend.order)

#heatmap showing % of signal in small size fraction 
ggplot(metab.heatmap.dat, aes(x = as.factor(Lat), y = reorder(Compound, -order), fill = Small_Frac_Perc)) +
  geom_tile() +
  scale_fill_viridis() +
  xlab("Latitude") +
  ylab("Compound")


















#heatmap showing % of signal in small size fraction 
ggplot(g2.all.dat, aes(x = as.factor(Lat), y = MF, fill = Sml_Frac_Perc)) +
  geom_tile() +
  scale_fill_viridis() +
  xlab("Latitude")

#bar chart showing distributions
ggplot(g2.all.dat, aes(x = as.factor(Lat), y = Mean.Area, fill = Size_Fraction)) +
  geom_col(position = "dodge") +
  facet_wrap(.~MF, scales = "free") +
  xlab("Latitude") +
  ylab("Peak Area")









