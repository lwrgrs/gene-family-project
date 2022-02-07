library(tidyverse)

setwd('~/Desktop/Calonectria_Projects/Gene_Families/Manuscript/R environments/')
load('annotation_objects_1.16.22')

# rename species names for plot
cog_classified_ogs$species_id <-
  recode(cog_classified_ogs$species_id,
          !!!c(Calonectriapseudonaviculata = "Calonectria pseudonaviculata",
            Calonectriahenricotiae = "Calonectria henricotiae",
            Calonectriamultiphialidica = "Calonectria multiphialidica",
            Calonectrianaviculata = "Calonectria naviculata",
            Calonectriapseudoreteaudii = "Calonectria pseudoreteaudii",
            Calonectrialeucothoes = "Calonectria leucothoes",
            Thelonectriarubi = "Thelonectria rubi",
            Neonectriahederae = "Neonectria hederae",
            Neonectriapunicea = "Neonectria punicea",
            Neonectriaditissima = "Neonectria ditissima",
            Corinectriafuckeliana = "Corinectria fuckeliana",
            Dactylonectriamacrodydima = "Dactylonectria macrodydima",
            Ilyonectriadestructans = "Ilyonectria destructans",
            Fusariumoxysporum4287 = "Fusarium oxysporum 4287",
            FusariumoxysporumFo47 = "Fusarium oxysporum Fo47",
            Fusariumfujikuroi = "Fusarium fujikuroi",
            Fusariumgraminearum = "Fusarium graminearum",
            Fusariumsolani = "Fusarium solani",
            Coccinonectriapachysandricola = "Coccinonectria pachysandricola",
            Pseudonectriabuxi = "Pseudonectria buxi",
            Stachybotryschlorohalonata = "Stachybotrys chlorohalonata",
            Stachybotryschartarum = "Stachybotrys chartarum"))

# set order of species on plot to reflect arrangement in phylogeny
level_order_dot <- c(
  "Calonectria pseudonaviculata",
  "Calonectria henricotiae",
  "Calonectria multiphialidica",
  "Calonectria naviculata",
  "Calonectria pseudoreteaudii",
  "Calonectria leucothoes",
  "Thelonectria rubi",
  "Neonectria hederae",
  "Neonectria punicea",
  "Neonectria ditissima",
  "Corinectria fuckeliana",
  "Dactylonectria macrodydima",
  "Ilyonectria destructans",
  "Fusarium oxysporum 4287",
  "Fusarium oxysporum Fo47",
  "Fusarium fujikuroi",
  "Fusarium graminearum",
  "Fusarium solani",
  "Coccinonectria pachysandricola",
  "Pseudonectria buxi",
  "Stachybotrys chlorohalonata",
  "Stachybotrys chartarum")

# ggplot of COG categories
dot_cog <- cog_classified_ogs %>%
  ungroup() %>%
  mutate(direction = case_when(
    gain_loss > 0 ~ "expanding",
    gain_loss < 0 ~ "contracting"
  )) %>%
  group_by(species_id, direction, cog_cat, OG) %>%
  summarise(prots_per_og = length(cog_cat)) %>%
  group_by(species_id, direction, cog_cat) %>%
  summarise(cog_cats_per_species = length(cog_cat)) %>%
  ggplot(aes(x = reorder(cog_cat, -cog_cats_per_species), 
             y = factor(species_id, levels = level_order_dot))) + 
  geom_point(aes(color = direction, size = cog_cats_per_species), alpha = 0.75) +
  scale_color_manual(values = c("#117733", "#332288"), 
                     labels = c("Rapidly Contracting", "Rapidly Expanding")) +
  scale_y_discrete(limits = rev) +
  scale_size(range = c(2, 10)) + # adjust range of point size
  labs(x = "COG Category", y = "Species",
       size = "Number of Rapidly Evolving\n Gene Families",
       color = "") +
  theme(axis.text.x = element_text(size = 9, color = "black"),
        axis.title.x = element_text(size = 11),
        axis.text.y = element_text(size = 9, 
                                   color = "black",
                                   face = 'italic'),
        axis.title.y = element_text(size = 11),
        legend.title = element_text(size = 10),
        legend.title.align = 0.5)

# write to file
setwd('~/Desktop/Calonectria_Projects/Gene_Families/Manuscript/Results/Figures/For_Submission/')
pdf('Figure_3_cbf_v2.pdf', width = 12, height = 7)

dot_cog

dev.off()