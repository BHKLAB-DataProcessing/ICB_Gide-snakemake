library(data.table)
library(readxl) 
library(stringr)
library(tximport)

args <- commandArgs(trailingOnly = TRUE)
work_dir <- args[1]
annot_dir <- args[2]

# CLIN_GIDE.txt
clin <- read_excel(file.path(work_dir, '1-s2.0-S1535610819300376-mmc2.xlsx'), sheet='Table S1. PD-1 Patient')
colnames(clin) <- clin[2, ]
clin <- clin[-c(1:2), ]
colnames(clin) <- str_replace_all(colnames(clin), '\\W', '.')
write.table(clin, file.path(work_dir, "CLIN_GIDE.txt"), sep="\t" , col.names=TRUE, row.names=FALSE)

# EXPR_gene_tpm.tsv, EXPR_gene_counts.tsv, EXPR_tx_tpm.tsv, EXPR_tx_counts.tsv
source('https://raw.githubusercontent.com/BHKLAB-Pachyderm/ICB_Common/main/code/process_kallisto_output.R')
load(file.path(annot_dir, "Gencode.v40.annotation.RData"))

dir.create(file.path(work_dir, 'rnaseq'))
zipfiles <- c('Gide_kallisto1.zip', 'Gide_kallisto2.zip', 'Gide_kallisto3.zip', 'Gide_kallisto4.zip', 'Gide_kallisto5.zip', 'Gide_kallisto6.zip')
for(zipfile in zipfiles){
  unzip(file.path(work_dir, zipfile), exdir=file.path(work_dir, 'rnaseq'))
}
unlink(file.path(work_dir, 'rnaseq', '__MACOSX'), recursive = TRUE)

process_kallisto_output(work_dir, tx2gene)