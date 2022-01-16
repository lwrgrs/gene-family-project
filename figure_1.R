library(tidyverse)
library(cowplot)

setwd('~/Desktop/Calonectria_Projects/Gene_Families/Manuscript/R environments/')
load('annotation_objects_1.16.22')

dat_means <-
  inner_join(means, dat, by = "species_id")

dat_means$species_id <-
  recode(dat_means$species_id,
          !!!c(Calonectriapseudonaviculata = "Calonectria pseudonaviculata",
            Calonectriahenricotiae = "Calonectria henricotiae",
            Calonectriamultiphialidica = "Calonectria multiphialidica",
            Calonectrianaviculata = "Calonectria naviculata",
            Calonectriapseudoreteaudii = "Calonectria pseudoreteaudii",
            Calonectrialeucothoes = "Calonectria leucothoes",
            Thelonectriarubi = "Thelonectria rubi",
            Aquanectriapenicillioides = "Aquanectria penicillioides",
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
            Pseudonectriafoliicola = "Pseudonectria foliicola",
            Stachybotryschlorohalonata = "Stachybotrys chlorohalonata",
            Stachybotryschartarum = "Stachybotrys chartarum"))

level_order <- c(
  "Calonectria pseudonaviculata",
  "Calonectria henricotiae",
  "Calonectria multiphialidica",
  "Calonectria naviculata",
  "Calonectria pseudoreteaudii",
  "Calonectria leucothoes",
  "Thelonectria rubi",
  "Aquanectria penicillioides",
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
  "Pseudonectria foliicola",
  "Stachybotrys chlorohalonata",
  "Stachybotrys chartarum")

dat_means <-
  dat_means %>%
  mutate(pct_rapid_ex = (n_expanding/(expns))*100,
         pct_rapid_co = (n_contracting/(contns))*100)

break_values <- pretty(dat_means_plot$n_gfs__)

dat_means_plot$direction <- reorder(dat_means_plot$direction, dat_means_plot$n_gfs__)

# expanding gene families
dat_means_plot_ex <- 
  dat_means %>%
  select(species_id, expns, n_expanding, pct_rapid_ex) %>%
  gather("direction_total", "n_total", c(expns, n_expanding))

dat_means_plot_ex$direction_total <- reorder(dat_means_plot_ex$direction_total, 
                                             dat_means_plot_ex$n_total)
dat_means_plot_ex$direction_total <- factor(dat_means_plot_ex$direction_total,
                                            levels(rev(dat_means_plot_ex$direction_total)))

# expanding gene families plot
e <- dat_means_plot_ex %>% 
  ggplot(aes(x=n_total, 
             y=factor(species_id, levels = level_order), 
             fill = direction_total)) +
  geom_bar(stat = "identity", 
           color = "black",
           position = "dodge") +
  scale_y_discrete(limits = rev) +
  scale_x_continuous(limits = c(-250, 12000)) +
  scale_fill_manual(values = c("#332288", "#88CCEE"),
                    labels = c("Rapidly Expanding", "Total Expanding")) +
  xlab("Number of Gene Families") +
  ylab("Species") +
  theme(legend.position = 'right',
        axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        legend.title = element_blank()) +
  geom_text(aes(label=paste(round(pct_rapid_ex, 2), "%"), x=pct_rapid_ex),
            hjust = 1.15)

# contracting gene families
dat_means_plot_co <- 
  dat_means %>%
  select(species_id, contns, n_contracting, pct_rapid_co) %>%
  gather("direction_total", "n_total", c(contns, n_contracting))

dat_means_plot_co$direction_total <- reorder(dat_means_plot_co$direction_total, 
                                             dat_means_plot_co$n_total)
dat_means_plot_co$direction_total <- factor(dat_means_plot_co$direction_total,
                                            levels(rev(dat_means_plot_co$direction_total)))

# contracting gene families plot
c <- dat_means_plot_co %>% 
  ggplot(aes(x=n_total, 
             y=factor(species_id, levels = level_order), 
             fill = direction_total)) +
  geom_bar(stat = "identity", 
           color = "black",
           position = "dodge") +
  scale_y_discrete(limits = rev) +
  scale_x_continuous(limits = c(-700, 12000)) +
  scale_fill_manual(values = c("#117733", "#999933"),
                    labels = c("Rapidly Contracting", "Total Contracting")) +
  xlab("Number of Gene Families") +
  ylab("Species") +
  theme(legend.position = 'right',
        axis.text.x = element_text(color = "black"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        legend.title = element_blank()) +
  geom_text(aes(label=paste(round(pct_rapid_co, 2), "%"), x=pct_rapid_co),
            hjust = 1.15)

# write each plot to a file (e and c to separate files)
setwd('~/Desktop/Calonectria_Projects/Gene_Families/Manuscript/Results/Figures/For_Submission/')

pdf('Figure_1_contracting_cbf.pdf', width = 13, height = 7)

c

dev.off()