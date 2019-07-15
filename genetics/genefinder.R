#!/usr/bin/env Rscript

# Libraries
library(tidyverse)
library(reshape2)

# Temporary pull in file
df <- read_tsv("phenotype-genotype-file.tab")

# Rename coronary dz to consolidate them
df$Trait <- gsub("Coronary Disease", "Coronary", df$Trait)
df$Trait <- gsub("Coronary Artery Disease", "Coronary", df$Trait)
df$Trait <- gsub("Acute Coronary Syndrome", "Coronary", df$Trait)
df$Trait <- gsub("Coronary Vasospasm", "Coronary", df$Trait)

# Rename Arrhythmias
df$Trait <- gsub("Death, Sudden, Cardiac", "Arrhythmia", df$Trait)
df$Trait <- gsub("Arrhythmias, Cardiac", "Arrhythmia", df$Trait)
df$Trait <- gsub("Heart Failure", "Arrhythmia", df$Trait)
df$Trait <- gsub("Cardiac Conduction Defect", "Arrhythmia", df$Trait)

# Rename Depression
df$Trait <- gsub("Depressive Disorder", "Depression", df$Trait)
df$Trait <- gsub("Depressive Disorder, Major", "Depression", df$Trait)
df$Trait <- gsub("Depression, Major", "Depression", df$Trait)
df$Trait <- gsub("Depression", "Depression", df$Trait)

# Rename Stress and Anxiety
df$Trait <- gsub("Anxiety Disorders", "Stress", df$Trait)
df$Trait <- gsub("Stress, Psychological", "Stress", df$Trait)
df$Trait <- gsub("Stress Disorders, Post-Traumatic", "Stress", df$Trait)
df$Trait <- gsub("Anxiety", "Stress", df$Trait)

# Important columns are selected and reorganized, stacking the Genes together
svar <- c("Trait", "SNP rs", "Gene", "Chromosome")
a <- df[names(df) %in% svar]
b <- df[c("Trait", "SNP rs", "Gene 2", "Chromosome")]
names(b)[3] <- "Gene"
df <- bind_rows(a,b) %>% unique()

# Only genes that are represented in more than 1 phenotype
df %<>%
  group_by(Gene) %>%
  filter(n() > 1)

# Widened the data, then finding only those with columns/traits of interest
dfw <- dcast(df[c("Trait", "Gene")], Gene ~ Trait) %>% as_tibble()

# Renumber all values to either 0 or 1
dfw[-1] -> tmp
tmp[tmp != 0] <- 1
dfw[-1] <- tmp

# look for rows with multiple positive findings
dfw$total <- rowSums(dfw[-1])
dfw <- dfw[dfw$total > 1, ]
dfw <- dfw[dfw$Depression != 0 | dfw$Stress != 0, ]
dfw <- dfw[dfw$Arrhythmia != 0 | dfw$Coronary != 0, ]

# Final table with SNPs and Chromosome position
tbl <- inner_join(dfw[-6], df[-1], by = "Gene") %>% unique() %>% group_by(Gene) %>% slice(1)
write.csv(tbl, "overlapping_genes.csv", row.names = FALSE)
