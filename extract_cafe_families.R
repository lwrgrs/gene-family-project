rm(list=ls())

library(readxl)
library(tidyverse)
library(seqinr)

# CREATE INDIVIDUAL FILES FOR EACH SPECIES CORRESPONDING TO ORTHOGROUPS AND RAPIDLY EVOLVING GENE FAMILIES

setwd("~/Desktop/Calonectria Projects/Gene Families/")

# read in orthofinder file
og <- read.csv("Orthogroups.csv", sep = "\t")
names <- c(colnames(og))

setwd("./Orthogroup Data/")

# make individual orthogroup files and corresponding protein sequences for each species
for (i in 2:25){
  write.csv(og[,c(1,i)], file = paste(names[i], "_og.csv", sep = ","))
}

setwd("~/Desktop/Calonectria_Projects/Gene_Families/")

# read in 
cafe_fams <- as.data.frame(read_xlsx("CAFE_fams_summary.xlsx", 3))
names <- c(colnames(cafe_fams))

setwd("~/Desktop/Calonectria Projects/Gene Families/Orthogroup Data/")

# create files that list rapidly evolving gene families for each species
for (i in 1:24){
  write.csv(cafe_fams[,i], file = paste(names[i], "se_og.csv"))
}

########
#######
######
#####
####
###
##
#

# CREATE FILES OF PROTEIN SEQUENCE IDS CORRESPONDING TO ORTHOGROUPS

rm(list=ls())

setwd("~/Desktop/Calonectria Projects/Gene Families/Orthogroup Data/")

# read in protein sequence IDs for each species
og_tru <- 
  as.data.frame(
    read.csv("Thelonectria_rubi_CBS113.12.proteins,_og.csv"))

# read in rapidly evolving gene family IDs for each species
cafe_fams_tru <- 
  as.data.frame(
    read.csv("Thelonectria rubi se_og.csv", 
             header = T, 
             row.names = 1))

# remove NAs
cafe_fams_tru <- 
  na.omit(cafe_fams_tru)

# create character vector of gene family IDs
se_og_tru <- 
  as.character(as.vector(cafe_fams_tru[,1]))

# select rows that match elements in the character vector of gene family IDs
selected_og_tru <- 
  og_tru[og_tru$X %in% se_og_tru,]

# concatenate protein seq IDs into a character vector
og_prot_seq_tru <- 
  paste(selected_og_tru[,3], collapse = ' ')

# look at the object, make sure there is nothing extra you need to do to
# clean it up
head(og_prot_seq_tru)

# split each element of the string into its own list element
ogprot_seq_tru <- 
  unlist(strsplit(og_prot_seq_tru, " "))

head(ogprot_seq_tru)

# remove empties
ogprotseq_tru <- 
  ogprot_seq_tru[ogprot_seq_tru != ""]

head(ogprotseq_tru)

# remove commas from each entry
ogps_tru <- 
  gsub(",$", "", ogprotseq_tru)

# final check before writing to file
ogps_tru

# write sequence names to file
write(ogps_tru, "~/Desktop/Calonectria Projects/Gene Families/Protein_sequences/OG_prot_seq_tru.txt")