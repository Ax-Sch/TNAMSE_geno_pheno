library(tidyverse)

# for debug
args=c("results/supporting2_anno/output_tsv.tsv", 
       "results/supporting2_anno/failed.txt", 
       "results/parse_table/all_monogenic_vars.tsv", 
       "../../")

args = commandArgs(trailingOnly=TRUE)
depth<-args[4]
converted_path<-paste0(depth, args[1]) #
failed_path<-paste0(depth, args[2])
unnested_path<-paste0(depth, args[3])

converted_vars<-read_csv2(converted_path)
failed_vars<-unlist(read_tsv(failed_path, col_names = FALSE, quote="'")) 
unnested_vars<-read_tsv(unnested_path) %>% 
  distinct()%>%
  left_join(converted_vars, by=c("HGVS_cDNA"="HGVS-Input")) %>%
  mutate(failed_conversion=HGVS_cDNA %in% failed_vars) %>%
  arrange(!failed_conversion)

write_tsv(file="unnested_vars_converted.tsv", x=unnested_vars)

patient_id_pos<-unnested_vars %>% 
  group_by(case_ID_paper)%>%
  summarise(chrom_positions=paste(ID, collapse=", "),
            HGVS=paste(HGVS_cDNA, collapse=", "),
                       HGVS_gDNA=paste(HGVS_gDNA, collapse=", "))

write_tsv(file="patient_id_pos.tsv", x=patient_id_pos)
