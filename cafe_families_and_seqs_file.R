rm(list=ls())

library(tidyverse)

setwd("~/Desktop/Calonectria_Projects/Gene_Families/Orthogroup Data/")

# read in protein sequence IDs corresponding to orthogroups for each species
che_og_prot <- 
  as.data.frame(
    read.csv("Calonectria_henricotiae_13.131.proteins,_og.csv"))

# read in list of rapidly evolving gene families for the species
che_cafe_og <- 
  as.data.frame(
    read.csv("Calonectria henricotiae se_og.csv", 
             header = T, 
             row.names = 1)) %>%
  na.omit(che_cafe_og)

# create vector of rapidly evolving gene family IDs for the species
che_og <- 
  as.character(as.vector(che_cafe_og[,1]))

# use vector to select orthogroups and associated protein sequences from che_og_prot data frame
che_og_prot_selected <- 
  che_og_prot[che_og_prot$X %in% che_og,]

# split string of protein seq IDs for each orthogroup at ","
che_og_prot_selected$Calonectria_henricotiae_13.131.proteins <-
  strsplit(che_og_prot_selected$Calonectria_henricotiae_13.131.proteins, ", ")

# unlist the protein seq IDs so they are separated into rows
che_og_prot_2 <-
  unnest(che_og_prot_selected, Calonectria_henricotiae_13.131.proteins) %>%
  as.data.frame(che_og_prot_2)

# write data frames to files to combine into one large file
write.csv(che_og_prot_2, "~/Desktop/Calonectria_Projects/Gene_Families/Orthogroup Data/REGF_sequences/che_regf_seq.csv")

# create master file of all rapidly evolving gene families and corresponding protein sequence IDs across species
all_regf <- read_csv("combined_regf_seq.csv",
                     col_names = T)

all_regf <- all_regf[complete.cases(all_regf),]
all_regf <- all_regf[,-c(1:2)]
colnames(all_regf) <- c("OG", "query")