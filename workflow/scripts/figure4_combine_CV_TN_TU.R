library(tidyverse)
library(readxl)

# for debugging 
args=c("results/figure4/ClinVar_table.tsv",
"results/figure4/TN_vars_count.tsv",
"results/figure4/Turro_variant_table.txt",
"resources/cancer_genes.txt",
"resources/genes_year_of_first_report.xlsx",
"resources/HGNC_symbols.txt",
"../../")

args = commandArgs(trailingOnly=TRUE)

depth<-args[7]
clinvar_path<-paste0(depth, args[1]) 
tnamse_path<-paste0(depth, args[2]) 
turro_path<-paste0(depth, args[3]) 
cancer_gene_path<-paste0(depth, args[4]) 
year_path<-paste0(depth, args[5]) 
hgnc_path<-paste0(depth, args[6]) 

clinvar<-read_tsv(clinvar_path)
tnamse<-read_tsv(tnamse_path)
turro<-read_tsv(turro_path)
cancer_gene<-read_tsv(cancer_gene_path)
year_tbl<-read_excel(year_path)
hgnc<-read_tsv(hgnc_path, col_names="GeneSymbol")$GeneSymbol

check_if_in_hgnc<-function(dataset, name_of_dataset){
  # find GeneSymbols not in hgnc and write to file
  not_in_hgnc<-dataset %>% filter (!GeneSymbol %in% hgnc)
  write_tsv(x=not_in_hgnc, file=paste0("not_in_HGNC_", name_of_dataset, ".tsv"))
  
  # give back a filtered list of variants in HGNC
  in_hgnc<-dataset %>% 
    filter (GeneSymbol %in% hgnc)%>%
    distinct(GeneSymbol, .keep_all = TRUE)
  return(in_hgnc)
}

TO_KEEP<-c("GeneSymbol", "number_of_variants")

# make sure all genes are in HGNC, just keep relevant columns
clinvar<-check_if_in_hgnc(clinvar, "ClinVar") %>% 
  select(any_of(TO_KEEP))
tnamse<-check_if_in_hgnc(tnamse, "TNamse")%>% 
  select(any_of(TO_KEEP))
turro<-check_if_in_hgnc(turro, "Turro")%>% 
  select(any_of(TO_KEEP))
year_tbl<-check_if_in_hgnc(year_tbl, "year_tbl")%>% 
  select(GeneSymbol, Year)
cancer_gene<-check_if_in_hgnc(cancer_gene, "cancer_gene")

# get all genes present in data sets
all_genes<-unique(c(
  clinvar$GeneSymbol,
  tnamse$GeneSymbol,
  turro$GeneSymbol
)) 

gene_counts_years<-tibble(GeneSymbol=all_genes)

# combine datasets
gene_counts_years<-gene_counts_years %>% 
  left_join(clinvar, by="GeneSymbol")%>%
  left_join(tnamse, by="GeneSymbol", suffix=c("","_tnamse"))%>%
  left_join(turro, by="GeneSymbol", suffix=c("_clinvar","_turro"))%>%
  left_join(year_tbl, by="GeneSymbol", suffix=c("","_year"))%>%
  mutate(cancer_gene = (GeneSymbol %in% unlist(cancer_gene,use.names = FALSE)))
  
head(gene_counts_years%>% arrange(-number_of_variants_clinvar))

# write to file
write_tsv(file="combined_gene_counts_years.tsv", gene_counts_years)
