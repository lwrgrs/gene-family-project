library(tidyverse)
library(ggtree)
library(tidytree)
library(treeio)
library(aplot)
library(scales)

setwd('~/Desktop/Calonectria_Projects/Gene_Families/Manuscript/R environments/')
load('annotation_objects_1.16.22')

setwd("~/Desktop/Calonectria_Projects/Gene_Families")

# read in phylogeny used in cafe analysis
tree <- read.tree("cafe_labeled_tree.tree") 

setwd("~/Desktop/Calonectria_Projects/Gene_Families/CAFE results/reports/")

dat_fig <- read_tsv("summary_run3_node.txt", col_names = T)
colnames(dat_fig) <- c("label", "expns", "contns", "regfs")

genus_tree <- c('Ilyonectria', 'Dactylonectria', 'Neonectria', 'Neonectria',
                'Neonectria', 'Corinectria', 'Calonectria', 'Calonectria',
                'Calonectria', 'Calonectria', 'Calonectria', 'Calonectria',
                'Thelonectria', 'Aquanectria', 'Fusarium', 'Fusarium',
                'Fusarium', 'Fusarium', 'Fusarium', 'Coccinonectria',
                'Pseudonectria', 'Pseudonectria', 'Stachybotrys', 'Stachybotrys')

species_tree <- c('destructans', 'macrodidyma', 'punicea', 'hederae',
                  'ditissima', 'fuckeliana', 'leucothoes', 'pseudoreteaudii',
                  'henricotiae', 'pseudonaviculata', 'multiphialidica', 'naviculata',
                  'rubi', 'penicillioides', 'oxysporum', 'oxysporum',
                  'fujikuroi', 'graminearum', 'solani', 'pachysandricola',
                  'foliicola', 'buxi', 'chartarum', 'chlorohalonata')

isolate_tree <- c(rep(NA, 14), "Fo47", "4287", rep(NA, 8))

# phylogeny with annotated branches
d <- data.frame(label = tree$tip.label, genus = genus_tree,
                species = species_tree, isolate = isolate_tree)

y <- full_join(tree, dat_fig, by  = 'label')
y <- as.treedata(y)

# additional cafe summary file with more tree annotations
setwd("~/Desktop/Calonectria_Projects/Gene_Families/CAFE results/reports/")
pub_dat <- read_tsv("summary_run3_pub.txt", col_names = T)

genus_pub <- c('Calonectria', 'Aquanectria', 'Calonectria', 'Stachybotrys',
               'Ilyonectria', 'Fusarium', 'Neonectria', 'Fusarium',
               'Dactylonectria', 'Calonectria', 'Coccinonectria', 'Pseudonectria',
               'Calonectria', 'Fusarium', 'Corinectria', 'Fusarium',
               'Neonectria', 'Stachybotrys', 'Fusarium', 'Pseudonectria', 
               'Neonectria', 'Thelonectria', 'Calonectria', 'Calonectria')

species_pub <- c('multiphialidica', 'penicillioides', 'henricotiae', 'chlorohalonata', 
                 'destructans', 'graminearum', 'hederae', 'oxysporum', 
                 'macrodidyma', 'leucothoes', 'pachysandricola', 'buxi',
                 'pseudonaviculata', 'solani', 'fuckeliana', 'oxysporum', 
                 'punicea', 'chartarum', 'fujikuroi', 'foliicola', 
                 'ditissima', 'rubi', 'naviculata', 'pseudoreteaudii')

isolate_pub <- c(rep(NA, 7), "Fo47", rep(NA, 7), "4287", rep(NA, 8))

p <- data.frame(expanded_fams = pub_dat$`Expanded fams`,
                contract_fams = pub_dat$`Contracted fams`,
                genus = genus_pub,
                species = species_pub,
                isolate = isolate_pub)

p$expanded_fams <- gsub(".*\\(", "", p$expanded_fams)
p$expanded_fams <- gsub("\\)", "", p$expanded_fams)
p$contract_fams <- gsub(".*\\(", "", p$contract_fams)
p$contract_fams <- gsub("\\)", "", p$contract_fams)

d <- inner_join(d, p, by = "species")
d <- d[,-c(7:8)]
d <- d[-c(16:17),]

colnames(d) <- c("label", "genus", "species",
                 "isolate",
                 "Rapidly Expanding", 
                 "Rapidly Contracting")

# mean gene gains and losses 
means <- og_seq_e_p %>%
  ungroup() %>%
  filter(gain_loss != 0) %>%
  mutate(direction = case_when(
    gain_loss > 0 ~ "expanding",
    gain_loss < 0 ~ "contracting"
  )) %>%
  group_by(species_id) %>%
  distinct(OG, .keep_all = T) %>%
  group_by(species_id, direction) %>%
  summarise(mean_gain_loss = round(mean(gain_loss), 2)) %>%
  spread(direction, mean_gain_loss, 2:3)

means$contracting[c(4,15,17,20)] <- 0
means$expanding[c(9,24)] <- 0
means$expanding[-c(9,24)] <- gsub("^", "\\+", means$expanding[-c(9,24)])

d <- d %>%
  arrange(label) %>%
  mutate(mean_loss = means$contracting,
         mean_gain = means$expanding) #%>%
unite("gain_loss", c(mean_gain, mean_loss), sep = " \\ ")

d_long <- d %>%
  gather("Rapidly Expanding", 
         "Rapidly Contracting", 
         key = direction, 
         value = total_change)

d_long$total_change <- as.numeric(d_long$total_change)

###
### assembling figure
###

# ggtree plot
t <- ggtree(tree, branch.length = 'none') %<+% d %<+% dat_fig + xlim(NA, 18) +
  geom_tiplab(aes(label = paste0('italic(', genus, 
                                 ')~bolditalic(', species, ')~', isolate)), 
              parse = TRUE,
              hjust = -0.05) +
  geom_tippoint(aes(size = regfs), color = "#999999", alpha = 1/2) +
  geom_nodepoint(aes(size = regfs), color = "#999999", alpha = 1/2) +
  theme(legend.position = "right",
        legend.title = element_text(size = 10),
        legend.title.align = 0.5) +
  labs(size = "Number of Rapidly Evolving\n Gene Families")

# geom_bar plot of rapidly expanding and contracting gene family counts
b1 <- ggplot(d_long, aes(fill = direction, 
                     x = total_change, 
                     y = label)) + 
  geom_bar(stat = 'identity',
           position = 'fill') + 
  scale_fill_manual(values = c("#117733", "#332288")) +
  scale_x_continuous(labels = percent) +
  theme(legend.position = 'right',
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        legend.title = element_blank(),
        axis.title.x = element_blank())

# join plots into one
b1 %>% insert_left(t, width = 2)

# write figure to file
setwd('~/Desktop/Calonectria_Projects/Gene_Families/Manuscript/Results/Figures/For_Submission')

pdf('Figure_2_cbf.pdf', width = 13, height = 7)

b1 %>% insert_left(t, width = 2)

dev.off()