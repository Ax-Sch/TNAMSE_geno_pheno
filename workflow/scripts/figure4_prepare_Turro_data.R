library(tidyverse)
library(readxl)

#development
args=c(
  "resources/41586_2020_2434_MOESM5_ESM.xlsx",
  "../../"
)
#setwd("../results/figure4/")

args = commandArgs(trailingOnly=TRUE)
depth<-args[2]
in_file<-paste0(depth, args[1]) #

snvs<-read_excel(path=in_file, sheet="SNV and indel list")
snvs_distinct<-snvs %>% 
  distinct(Gene, `Reference allele`, `Alternative allele`, `Position (GRCh37)`)

deletions<-read_excel(path=in_file, sheet="Large deletion list")

SVs<-read_excel(path=in_file, sheet="Complex structural variant list")

#GeneSymbol = na.omit(c(snvs$Gene, deletions$`Diagnostic-grade gene for the relevant domain`, SVs$`Diagnostic-grade gene for the relevant domain`))
## just use deletions that affect one gene:
deletions_one_gene_distinct<-deletions %>% filter(`Number of genes (all biotypes)` == 1)


GeneSymbol = na.omit(c(snvs_distinct$Gene, deletions_one_gene_distinct$`Diagnostic-grade gene for the relevant domain`))


gene_counts<-as.data.frame(table(GeneSymbol)) %>%
  rename(number_of_variants=Freq)

head(gene_counts)

write_tsv(file="Turro_variant_table.txt", gene_counts)