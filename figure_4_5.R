library(tidyverse)
library(aplot)

setwd('~/Desktop/Calonectria_Projects/Gene_Families/Manuscript/R environments/')
load('phi_figure.4.26.21.RData')

# rapidly evolving gene families with PHI annotated proteins
phi_gfs <-
  all_annot_bb %>% # phi annotations
  filter(grepl("^loss_of_pathogenicity$", `PHI Mutant Phenotype`) == T |
           grepl("^reduced_virulence$", `PHI Mutant Phenotype`) == T |
           grepl("effector", `PHI Mutant Phenotype`) == T) %>%
  distinct(`Rapidly Evolving Gene Families`)

# panel A of PHI figure (figure 4)
a_phi <- all_annot_bb %>%
  filter(grepl("^loss_of_pathogenicity$", `PHI Mutant Phenotype`) == T |
           grepl("^reduced_virulence$", `PHI Mutant Phenotype`) == T |
           grepl("effector", `PHI Mutant Phenotype`) == T) %>%
  mutate(`PHI Mutant Phenotype` = 
           fct_recode(`PHI Mutant Phenotype`,
                      "Loss of Pathogenicity" = "loss_of_pathogenicity",
                      "Reduced Virulence" = "reduced_virulence",
                      "Effector" = "effector_(plant_avirulence_determinant)",
                      "Effector" = "effector_(plant_avirulence_determinant)__unaffected_pathogenicity")) %>%
  group_by(`Rapidly Evolving Gene Families`, `PHI Mutant Phenotype`) %>%
  summarise(freq_phenotype = n()) %>%
  ggplot(aes(y = sort(`Rapidly Evolving Gene Families`, decreasing = F),
             x = freq_phenotype,
             fill = `PHI Mutant Phenotype`)) +
  geom_bar(stat="identity",
           position = "fill") +
  scale_x_continuous(labels = percent) +
  xlab("Percentage of PHI Annotated Sequences") +
  ylab("Rapidly Evolving Gene Families") +
  theme(legend.title = element_blank(),
        legend.position = "bottom",
        plot.title = element_text(hjust = 0.5),
        strip.text.y.left = element_text(angle = 0)) +
  scale_fill_manual(values = c("#44AA99", "#CC6677", "#DDCC77", "#CC6677"),
                    labels = c("E", "LOP",
                               "RV", 
                               "Not Annotated"))

# panel B of PHI figure (figure 4)
b_phi <- t_counts %>%
  filter(OG %in% phi_gfs$`Rapidly Evolving Gene Families`,
         grepl("^Calonectria", species) == T |
           grepl("Coccinonectriapachysandricola", species) == T |
           grepl("^Pseudonectria", species) == T) %>%
  mutate(species_2 = fct_recode(species,
                                "C. pseudonaviculata" = "Calonectriapseudonaviculata",
                                "C. henricotiae" = "Calonectriahenricotiae",
                                "C. multiphialidica" = "Calonectriamultiphialidica",
                                "C. naviculata" = "Calonectrianaviculata",
                                "C. leucothoes" = "Calonectrialeucothoes",
                                "C. pseudoreteaudii" = "Calonectriapseudoreteaudii",
                                "Co. pachysandricola" = "Coccinonectriapachysandricola",
                                "P. buxi" = "Pseudonectriabuxi",
                                "P. foliicola" = "Pseudonectriafoliicola"),
         direction = case_when(
           gain_loss > 0 ~ "Rapidly Expanding",
           gain_loss < 0 ~ "Rapidly Contracting",
           gain_loss == 0 ~ "Not Rapidly Evolving")) %>%
  ggplot(aes(y = OG, 
             x = factor(species_2, levels = level_order_og), 
             fill = direction)) +
  geom_tile(aes(width=0.7, height=0.7)) +
  theme(axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45,
                                   hjust = 1),
        legend.position = "bottom") +
  scale_fill_manual(values = c("#888888", "#117733", "#332288"))

# write figure 4 to file
setwd('~/Desktop/Calonectria_Projects/Gene_Families/Manuscript/Results/Figures/For_Submission/')
pdf('Figure_4_cbf.pdf',
    width = 11,
    height = 8)

plot_grid(a_phi, b_phi,
          labels = c("A", "B"),
          hjust = c("A" = -0.5, "B" = 0.4),
          label_size = 12)

dev.off()

# rapidly evolving gene families with predicted secreted protein sequences
sig_gfs <- all_annot_bb %>%
  filter(`Secreted Protein` == "Secreted" & `Effector Status` == "Effector" |
           `Secreted Protein` == "Secreted" & `Effector Status` == "Non-effector") %>%
  distinct(`Rapidly Evolving Gene Families`)

# panel A of figure 5
a_sig <- all_annot_bb %>%
  filter(`Secreted Protein` == "Secreted" & `Effector Status` == "Effector" |
           `Secreted Protein` == "Secreted" & `Effector Status` == "Non-effector") %>%
  group_by(`Rapidly Evolving Gene Families`, `Effector Status`) %>%
  summarise(freq_eff = n()) %>%
  ggplot(aes(y = sort(`Rapidly Evolving Gene Families`, decreasing = F),
             x = freq_eff,
             fill = `Effector Status`)) +
  geom_bar(stat="identity",
           position = "fill") +
  scale_x_continuous(labels = percent) +
  xlab("Percentage of Secreted Proteins Classified as Effectors") +
  ylab("Rapidly Evolving Gene Families") +
  theme(legend.title = element_blank(),
        legend.position = "bottom",
        plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values = c("#44AA99", "#DDCC77"))

# panel B of figure 5
b_sig <- t_counts %>%
  filter(OG %in% sig_gfs$`Rapidly Evolving Gene Families`,
         grepl("^Calonectria", species) == T |
           grepl("Coccinonectriapachysandricola", species) == T |
           grepl("^Pseudonectria", species) == T) %>%
  mutate(species_2 = fct_recode(species,
                                "C. pseudonaviculata" = "Calonectriapseudonaviculata",
                                "C. henricotiae" = "Calonectriahenricotiae",
                                "C. multiphialidica" = "Calonectriamultiphialidica",
                                "C. naviculata" = "Calonectrianaviculata",
                                "C. leucothoes" = "Calonectrialeucothoes",
                                "C. pseudoreteaudii" = "Calonectriapseudoreteaudii",
                                "Co. pachysandricola" = "Coccinonectriapachysandricola",
                                "P. buxi" = "Pseudonectriabuxi",
                                "P. foliicola" = "Pseudonectriafoliicola"),
         direction = case_when(
           gain_loss > 0 ~ "Rapidly Expanding",
           gain_loss < 0 ~ "Rapidly Contracting",
           gain_loss == 0 ~ "Not Rapidly Evolving")) %>%
  ggplot(aes(y = OG, 
             x = factor(species_2, levels = level_order_og), 
             fill = direction)) +
  geom_tile(aes(width=0.7, height=0.7)) +
  theme(axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45,
                                   hjust = 1),
        legend.position = "bottom") +
  scale_fill_manual(values = c("#888888", "#117733", "#332288"))

# write figure 5 to file
setwd('~/Desktop/Calonectria_Projects/Gene_Families/Manuscript/Results/Figures/For_Submission/')
pdf('Figure_5_cbf.pdf',
    width = 11,
    height = 8)

plot_grid(a_sig, b_sig,
          labels = c("A", "B"),
          hjust = c("A" = -0.5, "B" = 0.4),
          label_size = 12)

dev.off()