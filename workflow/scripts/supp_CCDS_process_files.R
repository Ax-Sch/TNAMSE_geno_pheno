# files downloaded from https://ftp.ncbi.nih.gov/pub/CCDS/current_human/
# 22 July 2021
library(tidyverse)

# for debugging
args=c("resources/CCDS2Sequence.current.txt",
       "resources/CCDS.current.txt",
       "resources/HGNC_gene_names_ids",
       "../../")

args = commandArgs(trailingOnly=TRUE)
depth<-args[4]
ressource_ccdsseq<-paste0(depth, args[1]) #
ressource_current<-paste0(depth, args[2]) 
hgnc_gene_names<-paste0(depth, args[3]) 


current_ccds_ids<-read_tsv(ressource_ccdsseq) %>%
  filter(current_member==1) %>% 
  distinct(`#ccds`) %>% 
  unlist()

ccds_original<-read_tsv(ressource_current, guess_max = 20000)

ccds_original$ccds_status %>% unique()

ccds_original<- ccds_original %>% filter(ccds_id %in% current_ccds_ids) 

ccds_original$ccds_status %>% unique()

ccds_reformated<-ccds_original %>%
  mutate(cds_locations=str_replace(cds_locations, "\\[", "")) %>%
  mutate(cds_locations=str_replace(cds_locations, "\\]", "")) %>%
  mutate(cds_locations=strsplit(cds_locations, ", ")) %>%
  unnest_longer(cds_locations)%>%
  separate(cds_locations, sep="-", into=c("cds_start","cds_end")) %>%
  mutate(cds_start=as.integer(cds_start), cds_end=as.integer(cds_end)) %>%
  mutate(cds_length=cds_end-cds_start) 

ccds_lengths<-
  ccds_reformated %>% 
  group_by(ccds_id,gene,gene_id,`#chromosome`,ccds_status) %>%
  summarise(ccds_total=sum(cds_length))%>%
  filter(ccds_status=="Public")%>%
  group_by(gene,gene_id,`#chromosome`,ccds_status)%>%
  summarise(max_ccds_gene=max(ccds_total))

bed_format<-ccds_reformated %>% 
  mutate(chrom=paste0("chr",`#chromosome`)) %>%
  select(`chrom`, cds_start, cds_end ) %>% 
  arrange(chrom,cds_start)%>%
  distinct()

write_tsv(bed_format, file="bedfile_total_cds.bed", col_names = FALSE)

HGNC_ids<-read_tsv(hgnc_gene_names)

ccds_lengths<-ccds_lengths %>% left_join(HGNC_ids, by=c("gene_id"="NCBI Gene ID(supplied by NCBI)")) %>%
  rename( "HGNC_name" =`Approved symbol`) %>% select(gene, gene_id,max_ccds_gene, HGNC_name)

write_tsv(ccds_lengths, file="ccds_lengths.tsv")
