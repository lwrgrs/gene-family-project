rm(list = ls())

library(tidyverse)
library(ggplot2)

# read in file containing total gene gains and losses for each species
setwd("~/Desktop/Calonectria_Projects/Gene_Families/CAFE results/reports/")

dat <- read_tsv("summary_run3_node.txt", col_names = T)
colnames(dat) <- c("label", "expns", "contns", "regfs")

dat <-
  dat %>%
  filter(grepl("[A-Z][a-z]*", label)) %>%
  mutate(species_id = gsub("<[0-9]*>$", "", label)) %>%
  select(species_id, expns, contns, regfs)

dat$n_expanding <-
  c(68, 1, 10, 5, 27, 1, 3, 86, 14, 21, 1, 1, 
    6, 49, 0, 169, 26, 17, 10, 1, 13, 0, 73, 14)

dat$n_contracting <-
  c(0, 3, 33, 13, 2, 9, 37, 20, 1, 4, 5, 0, 29,
    0, 11, 9, 5, 3, 46, 1, 0, 2, 1, 6)

# read in protein sequences corresponding to REGFs
setwd("~/Desktop/Calonectria_Projects/Gene_Families/Orthogroup Data/REGF_sequences/")

all_regf <- read_csv("combined_regf_seq.csv",
                     col_names = T)

positions <- which(all_regf$X.1 == "X.1")
pos <- c((5-1), (16-5), (122-16), (364-122), (1376-364), (2378-1376),
         (2389-2378), (2499-2389), (2693-2499), (2839-2693), (3017-2839),
         (3035-3017), (5011-3035), (5950-5011), (6469-5950), (6678-6469),
         (6765-6678), (6817-6765), (6989-6817), (6994-6989), (7003-6994),
         (7175-7003), (7244-7175), (7246-7244))
species_id <- all_regf[!complete.cases(all_regf),][4]

species <- c("Aquanectriapenicillioides", "Corinectriafuckeliana",
             "Calonectriahenricotiae", "Calonectrialeucothoes",
             "Calonectriamultiphialidica", "Calonectrianaviculata",
             "Coccinonectriapachysandricola", "Calonectriapseudonaviculata",
             "Calonectriapseudoreteaudii", "Dactylonectriamacrodydima",
             "Fusariumfujikuroi", "Fusariumgraminearum",
             "Fusariumoxysporum4287", "FusariumoxysporumFo47", 
             "Fusariumsolani", "Ilyonectriadestructans",
             "Neonectriaditissima", "Neonectriahederae",
             "Neonectriapunicea", "Pseudonectriabuxi",
             "Pseudonectriafoliicola", "Stachybotryschartarum",
             "Stachybotryschlorohalonata", "Thelonectriarubi")

v <- c()
for (i in 1:24){
  values <- rep(species[i], pos[i])
  v <- c(v, values)
}

all_regf <- all_regf[-5,]
all_regf$species <- unlist(v)
all_regf <- all_regf[complete.cases(all_regf),]
all_regf <- all_regf[,-c(1:2)]
colnames(all_regf) <- c("OG", "query", "species")

rm(species_id, v, values, i, pos, positions) # clear out environment

# read in count summary for cafe families
setwd("~/Desktop/Calonectria_Projects/Gene_Families/CAFE results/Summaries/")

counts <- read.csv("gene_family_counts.csv", header = T, row.names = 1)

t_counts <- as.data.frame(t(counts)) %>%
  tibble::rownames_to_column("OG") %>%
  gather("species", "gain_loss", 2:25)

# read in and format eggnog annotations
setwd("~/Desktop/Calonectria_Projects/Gene_Families/eggNOG/")

c_egg_annot <- read_tsv("combined_eggnog_annot.tsv",
                        col_names = F,
                        skip = 4,
                        comment = "#")

c_egg_annot$X21 <- as.factor(c_egg_annot$X21)

c_eggnog <- data.frame(query = c_egg_annot$X1, 
                       seed_evalue = c_egg_annot$X3, 
                       cog_cat = c_egg_annot$X21)

# read in and format pfam annotations
setwd("~/Desktop/Calonectria_Projects/Gene_Families/OG Annotations/Annotated_files/Pfam_annot/")

c_pfam_annot <- read_csv("combined_pfam_annotations.csv",
                         col_names = T)

c_pfam_annot <- c_pfam_annot[complete.cases(c_pfam_annot),]
c_pfam_annot <- c_pfam_annot[,-c(1,9)]
colnames(c_pfam_annot) <- c("pfam_target",
                            "pfam_accession",
                            "query",
                            "pfam_evalue",
                            "hmm_score",
                            "hmm_bias",
                            "target_description",
                            "OG")

c_pfam_filter <- 
  c_pfam_annot %>%
  group_by(query) %>%
  filter(pfam_evalue == min(pfam_evalue))

####
#### match annotations to protein sequences and rapidly evolving gene families in a single data frame
####

# join protein sequence IDs with orthogroups and gain/loss counts (t_counts)
og_seq <-
  right_join(all_regf, t_counts, by = c("OG", "species"))

# join eggnog annotations
og_seq_e <-
  right_join(c_eggnog, og_seq, by = "query")

# join pfam annotations
og_seq_e_p <-
  right_join(c_pfam_filter, og_seq_e, by = "query")

# remove duplicated columns
og_seq_e_p <- og_seq_e_p[,-c(5:8)]
colnames(og_seq_e_p) <- c("pfam_target", "pfam_accession",
                          "query", "pfam_evalue", "eggnog_evalue", "cog_cat", "OG",
                          "species_id", "gain_loss")

og_seq_e_p[og_seq_e_p$cog_cat == "<NA>",] <- NA

# join protein counts
og_seq_e_p <-
  right_join(dat, og_seq_e_p, by = "species_id")

# filter sequences with NA pfam and cog annotations
og_seq_ep_pfam_filter <-
  og_seq_e_p %>%
  ungroup() %>%
  distinct(query, .keep_all = TRUE) %>%
  filter(is.na(pfam_target) == F)

# data frame of annotated cafe families and corresponding protein sequences
cog_classified_ogs <-
og_seq_e_p %>%
  ungroup() %>%
  distinct(query, .keep_all = TRUE) %>%
  filter(is.na(cog_cat) == F) %>%
  group_by(OG, cog_cat) %>%
  summarise(max_cog_cat = length(cog_cat)) %>%
  filter(max_cog_cat == max(max_cog_cat)) %>%
filter(OG != "OG0000156" | cog_cat != "S",
       OG != "OG0000612" | cog_cat != "S",
       OG != "OG0000648" | cog_cat != "I",
       OG != "OG0000649" | cog_cat != "M",
       OG != "OG0001150" | cog_cat != "D") %>%
  select(OG, cog_cat) %>%
  inner_join(og_seq_ep_pfam_filter, by = c("OG", "cog_cat"))