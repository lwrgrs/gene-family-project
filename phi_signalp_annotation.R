rm(list=ls())

library(tidyverse)
library(Biostrings)
library(readxl)

setwd("~/Desktop/Calonectria_Projects/Gene_Families/PHI/blastp_results/")

# read in blastp results for each species
che_phi <- read_tsv("phi_blast_che.txt", col_names = F)
cps_phi <- read_tsv("phi_blast_cps.txt", col_names = F)
cle_phi <- read_tsv("phi_blast_cle.txt", col_names = F)
cmu_phi <- read_tsv("phi_blast_cmu.txt", col_names = F)
cna_phi <- read_tsv("phi_blast_cna.txt", col_names = F)
cpseudor_phi <- read_tsv("phi_blast_cpseudor.txt", col_names = F)
cpa_phi <- read_tsv("phi_blast_cpa.txt", col_names = F)
pbu_phi <- read_tsv("phi_blast_pbu.txt", col_names = F)
pfo_phi <- read_tsv("phi_blast_pfo.txt", col_names = F)

phi_hits <- rbind(che_phi, cps_phi, cle_phi, cmu_phi, cna_phi,
                  cpseudor_phi, cpa_phi, pbu_phi, pfo_phi)

blast_colnames <- c("query", "subject", "pct_matches", "length",
                    "mismatch", "gaps", "qstart", "qend", "sstart", "send",
                    "evalue", "bitscore")

colnames(phi_hits) <- blast_colnames

# remove unneeded objects
rm(che_phi, cps_phi, cle_phi, cmu_phi, cna_phi, cpseudor_phi, cpa_phi, pbu_phi, pfo_phi)

phi_hits_filtered <-
  phi_hits %>%
  group_by(query) %>%
  filter(evalue == min(evalue)) %>%
  inner_join(all_regf, by = "query") %>%
  separate(subject, into = c("phi_id",
                             "phi_accession",
                             "gene",
                             "phi_number",
                             "pathogen_species",
                             "phenotype"),
           sep = "#")

# signalp and effectorp annotations
setwd("~/Desktop/Calonectria_Projects/Gene_Families/SignalP/")

# read in signalp predictions
combined_sp <- read_tsv("combined_signalp_hits.txt", 
                        col_names = F,
                        comment = "#")
colnames(combined_sp) <- c("query", "prediction", "sp", "other", "cs_position")

og_sigp <-
  inner_join(combined_sp, all_regf, by = "query")

# read in effectorp predictions
effp <-
  read_xlsx("~/Desktop/Calonectria_Projects/Gene_Families/EffectorP/effectorp_results.xlsx",
                    col_names = T)
colnames(effp) <- c("query", "effector", "probability")

og_sigp_effp <-
  right_join(effp, og_sigp, by = "query")