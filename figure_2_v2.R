library(tidyverse)
library(ggtree)
library(tidytree)
library(treeio)
library(aplot)
library(scales)
library(patchwork)

setwd('~/Desktop/Calonectria_Projects/Gene_Families/Manuscript/R environments/')
load('annotation_objects_1.16.22')

# read in phylogeny used in cafe analysis and additional tree annotations
setwd("~/Desktop/Calonectria_Projects/Gene_Families/CAFE_input_data/")
tree <- read.tree("ultrameric_tree_V1_cafe") 

# total rapidly evolving gene family annotations
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

# matching annotations to tree labels
d <- data.frame(label = tree$tip.label, genus = genus_tree,
                species = species_tree, isolate = isolate_tree)

# number of expanding vs. contracting gene family annotations
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

d <- d %>%
  arrange(label) %>%
  mutate(mean_loss = means$contracting,
         mean_gain = means$expanding)

d_long <- d %>%
  gather("Rapidly Expanding", 
         "Rapidly Contracting", 
         key = direction, 
         value = total_change)

d_long$total_change <- as.numeric(d_long$total_change)

# annotated ggtree
t <- ggtree(tree) %<+% d %<+% dat_fig +
  geom_tiplab(aes(label = paste0('italic(', genus, 
                                 ')~italic(', species, ')~', isolate)), 
              parse = TRUE,
              hjust = -0.05) + 
  geom_tippoint(aes(size = regfs), color = "#999999", alpha = 1/2) +
  geom_nodepoint(aes(size = regfs), color = "#999999", alpha = 1/2) +
  theme(legend.position = "right",
        legend.title = element_text(size = 10),
        legend.title.align = 0.5) +
  labs(size = "Number of Rapidly Evolving \nGene Families") +
  coord_cartesian(clip = 'off') +
  theme_tree2(plot.margin = margin(6, 190, 6, 6),
              legend.title.align = 0.5,
              legend.title = element_text(size = 10.5)) + 
  scale_x_continuous(labels = abs)

# bar graph of rapidly evolving gene family percentages
order_sp <- c("Stachybotryschartarum<44>", "Stachybotryschlorohalonata<46>",
              "Pseudonectriafoliicola<40>", "Pseudonectriabuxi<42>",
              "Coccinonectriapachysandricola<38>", "Fusariumsolani<36>",
              "Fusariumgraminearum<34>", "Fusariumfujikuroi<32>",
              "FusariumoxysporumFo47<28>", "Fusariumoxysporum4287<30>",
              "Ilyonectriadestructans<0>", "Dactylonectriamacrodydima<2>",
              "Corinectriafuckeliana<10>", "Neonectriaditissima<8>",
              "Neonectriapunicea<4>", "Neonectriahederae<6>",
              "Aquanectriapenicillioides<26>", "Thelonectriarubi<24>",
              "Calonectrialeucothoes<12>", "Calonectriapseudoreteaudii<14>",
              "Calonectrianaviculata<22>", "Calonectriamultiphialidica<20>",
              "Calonectriahenricotiae<16>", "Calonectriapseudonaviculata<18>")

b1 <- ggplot(d_long, aes(fill = direction, 
                         x = total_change, 
                         y = factor(label, levels = order_sp))) + 
  geom_bar(stat = 'identity',
           position = 'fill') + 
  scale_fill_manual(values = c("#117733", "#332288")) +
  scale_x_continuous(labels = percent) +
  theme(legend.position = 'right',
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 10.5),
        axis.title.x = element_blank())

# different ways to plot both
b1 %>% insert_left(revts(t), width = 2)

revts(t) + b1 + plot_layout(guides = 'collect')

# write to file
setwd('~/Desktop/Calonectria_Projects/Gene_Families/Manuscript/Results/Figures/For_Submission')

pdf('Figure_2_cbf_v2.pdf', width = 13, height = 10)

revts(t) + b1 + plot_layout(guides = 'collect')

dev.off()