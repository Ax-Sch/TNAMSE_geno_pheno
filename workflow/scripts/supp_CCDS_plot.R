library(tidyverse)
library(readxl)
library(ggpubr)
# for debug
args = c("results/parse_table/diagnostic_vars.tsv", 
         "results/supporting_ccds/ccds_lengths.tsv",
         "resources/mimTitles.txt", 
         "resources/genemap2_15_07_2021.txt", 
         "../../")

args = commandArgs(trailingOnly=TRUE)
depth<-args[5]
variants_TN_file<-paste0(depth, args[1]) #
ccds_lengths_file<-paste0(depth, args[2]) 
mim_titles_file<-paste0(depth, args[3]) # "resources/mimTitles.txt"
genemap2_file<-paste0(depth, args[4])

variants_TN<-read_tsv(variants_TN_file) %>%
  rename(path_classes_split=ACMG_class,
         Disease_gene_split=gene)

#candidate_genes<-read_xlsx(candidates_file)%>%
#  filter(!is.na(Gene_name))

#unique(candidate_genes$Classification)

DGG_or_high<-variants_TN%>%
  filter(candidate)

n_cases<-length(unique(DGG_or_high$case_ID_paper))
n_genes<-length(unique(DGG_or_high$Disease_gene_split))

print(paste0("Number of cases: ", n_cases))
print(paste0("Number of genes: ", n_genes))

DGG_or_high_gene_name<-unique(DGG_or_high$Disease_gene_split)

ccds_lengths<-read_tsv(ccds_lengths_file) %>% 
  filter(!is.na(HGNC_name))

# restrict on number sign OMIM entries
mim_titles<-read_tsv(mim_titles_file, skip=2) %>%
  filter(`# Prefix` == "Number Sign") %>% 
  distinct(`MIM Number`) %>%
  unlist()

genemap2<-read_tsv(genemap2_file, skip=3)


TN_diagnostic_genes <- variants_TN %>% 
  filter(!Disease_gene_split %in% DGG_or_high_gene_name)%>%
  distinct(Disease_gene_split) %>%  
  left_join(ccds_lengths, by=c("Disease_gene_split"="HGNC_name")) %>%
  mutate(case="TN_Diagnostic") %>% 
  dplyr::select(case, max_ccds_gene)

TN_novel_genes <- DGG_or_high %>% 
  distinct(Disease_gene_split) %>%  
  left_join(ccds_lengths, by=c("Disease_gene_split"="HGNC_name")) %>%
  mutate(case="TN_Research") %>% 
  dplyr::select(case, max_ccds_gene)

genemap2<-genemap2 %>%
  mutate(Phenotypes=strsplit(Phenotypes, ", ")) %>%
  unnest_longer(Phenotypes)%>%
  mutate(Phenotypes=str_extract(Phenotypes, "[0-9][0-9][0-9][0-9][0-9][0-9]")) %>%
  filter(!is.na(Phenotypes))

genemap_number<-genemap2 %>%
  filter(Phenotypes %in% mim_titles)%>%
  left_join(ccds_lengths, by=c("Entrez Gene ID"="gene_id")) %>%
  filter(!is.na(gene)) %>% 
  distinct(gene,max_ccds_gene)

#write_tsv(genemap2, file="result_for_manual_check_in_Excel.tsv")

####### do manual check !!!!

OMIM_genes<-genemap_number %>%
  mutate(case="OMIM") %>% 
  dplyr::select(case, max_ccds_gene)

all_genes<-ccds_lengths%>%
  mutate(case="All_genes") %>% 
  dplyr::select(case, max_ccds_gene)

joined_cases<-rbind(TN_diagnostic_genes,
                    TN_novel_genes,
                    OMIM_genes,
                    all_genes
) %>%
  filter(!is.na(max_ccds_gene))

print(table(joined_cases$case))

breaks<-c(10,100,300, 1000,3000,10000,30000, 100000)
p1<-ggplot(joined_cases,
           aes(x=case,y=max_ccds_gene)) +
  geom_jitter(width=0.25, color="dodgerblue2")+
  geom_boxplot(outlier.shape = NA, fill=NA, color="black")+
  theme_classic()+
  scale_y_log10(breaks = breaks, labels = breaks)+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.95, hjust=1))
  

print(p1)

ggsave("Gene_lengths_log.pdf", p1, width=5, height=5)

pairwise.t.test(joined_cases$max_ccds_gene, joined_cases$case, p.adjust.method="none")
pairwise.t.test(joined_cases$max_ccds_gene, joined_cases$case, p.adjust.method="bonferroni")

