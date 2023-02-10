library(tidyverse)

# for development
# setwd("../results/figure4/")
args=c(
  "results/figure4/combined_gene_counts_years.tsv",
  "../../"
)

args = commandArgs(trailingOnly=TRUE)
depth<-args[2]
all_w_clinvar_path<-paste0(depth, args[1]) #

all_w_clinvar<-read_tsv(all_w_clinvar_path)

# Prepare the data
all_w_clinvar_no_tumor<-all_w_clinvar %>%
  mutate(number_of_variants_clinvar=replace_na(data=number_of_variants_clinvar, replace=0))%>%
  filter(cancer_gene==FALSE)  %>%
  arrange(-number_of_variants_clinvar)%>% 
  mutate(cumsum_clinvar=cumsum(number_of_variants_clinvar))%>% 
  mutate(cumsum_TNAMSE=cumsum(number_of_variants_tnamse))%>%
  mutate(cumsum_Turro=cumsum(number_of_variants_turro))

total_turro<-sum(all_w_clinvar_no_tumor$number_of_variants_turro, na.rm = TRUE)
total_TNAMSE<-sum(all_w_clinvar_no_tumor$number_of_variants_tnamse, na.rm = TRUE)
total_clinvar<-sum(all_w_clinvar_no_tumor$number_of_variants_clinvar, na.rm=TRUE)

max_TNAMSE<-max(all_w_clinvar_no_tumor$number_of_variants_tnamse, na.rm = TRUE)
max_Turro<-max(all_w_clinvar_no_tumor$number_of_variants_turro, na.rm = TRUE)
max_clinvar<-max(all_w_clinvar_no_tumor$number_of_variants_clinvar, na.rm=TRUE)

max_cumsum_TNAMSE<-max(all_w_clinvar_no_tumor$cumsum_TNAMSE[!is.na(all_w_clinvar_no_tumor$cumsum_TNAMSE)])
max_cumsum_clinvar<-max(all_w_clinvar_no_tumor$cumsum_clinvar, na.rm=TRUE)


quart_clinv<-total_clinvar/4
print("total clinvar")
print(total_clinvar)

print("Clinvar quarter size")
print(quart_clinv)

all_w_clinvar_no_tumor<-all_w_clinvar_no_tumor %>% 
  mutate(clinvar_quarter=cut(all_w_clinvar_no_tumor$cumsum_clinvar, 
                             breaks=c(0, quart_clinv*1, quart_clinv*2, quart_clinv*3, quart_clinv*4), 
                             labels=c("first","second","third","fourth"))) %>%
  mutate(year_range=cut(Year, 
                        breaks=c(1900, 2000, 2010, 2022), 
                        labels=c("-2000","2001-2010","2011-2022")))%>%
  mutate(rel_tnamse=number_of_variants_tnamse/total_TNAMSE)%>%
  mutate(rel_turro=number_of_variants_turro/total_turro)%>%
  mutate(rel_clinvar=number_of_variants_clinvar/total_clinvar)


cv_quarters<-all_w_clinvar_no_tumor%>% 
  group_by(clinvar_quarter)%>%
  summarise(max_cv=max(cumsum_clinvar))%>%
  distinct(clinvar_quarter, max_cv)


print("genes in quarters")
print(table(all_w_clinvar_no_tumor$clinvar_quarter))
print("cumulative vars relative to quarters")
print(cv_quarters)

all_w_clinvar_no_tumor_sorted_year<-all_w_clinvar_no_tumor %>% 
  arrange(str_sort(year_range, numeric=TRUE))%>%
  mutate(year_range=factor(year_range, levels=unique(year_range)))

all_w_clinvar_no_tumor_sorted_clinvar<-all_w_clinvar_no_tumor %>% 
  arrange(desc(clinvar_quarter))%>%
  mutate(clinvar_quarter=factor(clinvar_quarter, levels=unique(clinvar_quarter)))%>%
  filter(!is.na(clinvar_quarter))


# 2D plots with ClinVar as Reference
freq_clinvar_vs_tnamse<-ggplot(all_w_clinvar_no_tumor)+
  geom_line(aes(x=cumsum_clinvar,y=number_of_variants_clinvar/max_clinvar), size=0.5)+
  geom_point(aes(x=cumsum_clinvar,y=number_of_variants_clinvar/max_clinvar), size=1)+
  geom_jitter(data=all_w_clinvar_no_tumor[!is.na(all_w_clinvar_no_tumor$number_of_variants_tnamse),], size=1.0, alpha=0.7, width=0.05,height=0.05, aes(x=cumsum_clinvar,y=number_of_variants_tnamse/max_TNAMSE, color=year_range))+
  theme_bw()
print(freq_clinvar_vs_tnamse)
ggsave(freq_clinvar_vs_tnamse, filename = "freq_clinvar_vs_tnamse.pdf", width=4.5, height=3.0)

freq_clinvar_vs_turro<-ggplot(all_w_clinvar_no_tumor)+
  geom_line(aes(x=cumsum_clinvar,y=number_of_variants_clinvar/max_clinvar), size=0.5)+
  geom_point(aes(x=cumsum_clinvar,y=number_of_variants_clinvar/max_clinvar), size=1)+
  geom_jitter(data=all_w_clinvar_no_tumor[!is.na(all_w_clinvar_no_tumor$number_of_variants_turro),], size=1.0, alpha=0.7, width=0.01,height=0.01, aes(x=cumsum_clinvar,y=number_of_variants_turro/max_Turro, color=year_range))+
  theme_bw()
print(freq_clinvar_vs_turro)
ggsave(freq_clinvar_vs_turro, filename = "freq_clinvar_vs_turro.pdf", width=4.5, height=3.0)





# Bar plots numbers
plot_number<-ggplot() + 
  geom_bar(data=all_w_clinvar_no_tumor %>% arrange(-number_of_variants_turro), aes(x="Turro et al.",y=number_of_variants_turro), color="grey40",position="stack", fill="white",stat="identity", width=0.7) + 
  geom_bar(data=all_w_clinvar_no_tumor %>% arrange(-number_of_variants_tnamse), aes(x="TNAMSE",y=number_of_variants_tnamse), color="grey40",position="stack", fill="white",stat="identity", width=0.7) + 
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45))
plot_number
ggsave(file="plot_numbers.pdf",plot_number, width=1.8, height=3.2)

# Bar plots years
plot_year_range_abs<-ggplot() + 
  geom_bar(data=all_w_clinvar_no_tumor_sorted_year %>% arrange(-number_of_variants_turro), aes(x="Turro et al.",y=number_of_variants_turro, color=year_range),position="stack", fill="white",stat="identity") + 
  geom_bar(data=all_w_clinvar_no_tumor_sorted_year %>% arrange(-number_of_variants_tnamse), aes(x="TNAMSE",y=number_of_variants_tnamse, color=year_range),position="stack", fill="white",stat="identity") + 
  theme_minimal()
plot_year_range_abs
ggsave(file="plot_year_range_abs.pdf",plot_year_range_abs, width=3, height=5)

plot_year_range_rel<-ggplot() + 
  geom_bar(data=all_w_clinvar_no_tumor_sorted_year %>% arrange(-number_of_variants_turro), aes(x="Turro et al.",y=rel_turro, color=year_range),position="stack", fill="white",stat="identity") + 
  geom_bar(data=all_w_clinvar_no_tumor_sorted_year %>% arrange(-number_of_variants_tnamse), aes(x="TNAMSE",y=rel_tnamse, color=year_range),position="stack", fill="white",stat="identity") + 
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45))
plot_year_range_rel
ggsave(file="plot_year_range_rel.pdf",plot_year_range_rel, width=3, height=5)

plot_year_turro<-ggplot(data=all_w_clinvar_no_tumor_sorted_year %>% arrange(-number_of_variants_turro) ) + #%>% filter(!is.na(year_range))
  geom_bar( aes(x=year_range,y=number_of_variants_turro, color=year_range),position="stack", fill="white",stat="identity", width=0.7) + 
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, hjust=0.95))
plot_year_turro
ggsave(file="plot_year_turro.pdf",plot_year_turro, width=4.0, height=3.2)

plot_year_tn<-ggplot(data=all_w_clinvar_no_tumor_sorted_year %>% arrange(-number_of_variants_tnamse)) + # %>% filter(!is.na(year_range))
  geom_bar( aes(x=year_range,y=number_of_variants_tnamse, color=year_range),position="stack", fill="white",stat="identity", width=0.7) + 
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, hjust=0.95))
plot_year_tn
ggsave(file="plot_year_tn.pdf",plot_year_tn, width=4.0, height=3.2)




# Bar plots ClinVar quartiles
ggplot() + 
  geom_bar(data=all_w_clinvar_no_tumor_sorted_clinvar %>% arrange(-number_of_variants_turro) %>% filter(!is.na(clinvar_quarter)), aes(x="Turro et al.",y=number_of_variants_turro, color=clinvar_quarter),position="stack", fill="white",stat="identity") + 
  geom_bar(data=all_w_clinvar_no_tumor_sorted_clinvar %>% arrange(-number_of_variants_tnamse) %>% filter(!is.na(clinvar_quarter)), aes(x="TNAMSE",y=number_of_variants_tnamse, color=clinvar_quarter),position="stack", fill="white",stat="identity")+
  theme_minimal()

plot_clinvar<-ggplot() + 
  geom_bar(data=all_w_clinvar_no_tumor_sorted_clinvar %>% arrange(-number_of_variants_turro), aes(x="Turro et al.",y=rel_turro, color=clinvar_quarter),position="stack", fill="white",stat="identity") + 
  geom_bar(data=all_w_clinvar_no_tumor_sorted_clinvar %>% arrange(-number_of_variants_tnamse), aes(x="TNAMSE",y=rel_tnamse, color=clinvar_quarter),position="stack", fill="white",stat="identity") + 
  geom_bar(data=all_w_clinvar_no_tumor_sorted_clinvar %>% arrange(-number_of_variants_clinvar), aes(x="ClinVar",y=rel_clinvar, color=clinvar_quarter),position="stack", fill="white",stat="identity") + 
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45))

plot_clinvar


plot_clinvar<-ggplot() + 
  geom_bar(data=all_w_clinvar_no_tumor_sorted_clinvar %>% arrange(-number_of_variants_tnamse), aes(x="TNAMSE",y=rel_tnamse, color=clinvar_quarter, fill=clinvar_quarter),position="stack", stat="identity") + 
  geom_bar(data=all_w_clinvar_no_tumor_sorted_clinvar %>% arrange(-number_of_variants_clinvar), aes(x="ClinVar",y=rel_clinvar, color=clinvar_quarter, fill=clinvar_quarter),position="stack",stat="identity") + 
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45))


ggsave(file="plot_clinvar_quartile_stacked.pdf",plot_clinvar, width=3, height=5)


plot_year_range<-ggplot() + 
  geom_bar(data=all_w_clinvar_no_tumor_sorted_clinvar %>% arrange(-number_of_variants_tnamse), aes(x=year_range,y=rel_tnamse, fill=year_range, color=year_range),position="stack", stat="identity") + 
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45))+
  facet_grid(vars(clinvar_quarter), scales="free")


ggsave(file="plot_year_range_binned_by_clinvar.pdf",plot_year_range, width=3, height=6)



ggsave(file="plot_clinvar.pdf",plot_clinvar, width=3, height=5)




all_w_clinvar_no_tumor_sorted_clinvar<-all_w_clinvar_no_tumor %>% 
  arrange((clinvar_quarter))%>%
  mutate(clinvar_quarter=factor(clinvar_quarter, levels=unique(clinvar_quarter)))

plot_clinvar_clinvar<-ggplot() + 
  geom_bar(data=all_w_clinvar_no_tumor_sorted_clinvar %>% arrange(-number_of_variants_clinvar), aes(x=clinvar_quarter,y=number_of_variants_clinvar, color=clinvar_quarter),position="stack", fill="white",stat="identity") + 
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45))

ggsave(file="plot_clinvar_clinvar.pdf",plot_clinvar_clinvar, width=4.5, height=7)


plot_clinvar_turro<-ggplot() + 
  geom_bar(data=all_w_clinvar_no_tumor_sorted_clinvar %>% arrange(-number_of_variants_turro), aes(x=clinvar_quarter,y=number_of_variants_turro, color=clinvar_quarter),position="stack", fill="white",stat="identity", width=0.7) + 
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, hjust=0.95))

ggsave(file="plot_clinvar_turro.pdf",plot_clinvar_turro, width=4.0, height=3.2)

plot_clinvar_tn<-ggplot() + 
  geom_bar(data=all_w_clinvar_no_tumor_sorted_clinvar %>% arrange(-number_of_variants_tnamse), aes(x=clinvar_quarter,y=number_of_variants_tnamse, color=clinvar_quarter),position="stack", fill="white",stat="identity", width=0.7) + 
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, hjust=0.95))

ggsave(file="plot_clinvar_tn.pdf",plot_clinvar_tn, width=4.0, height=3.2)

#### Kolmogorov-Smirnoff-Test
turro <- rep(unique(all_w_clinvar_no_tumor$Year),times=sapply(unique(all_w_clinvar$Year),function(x){sum(all_w_clinvar_no_tumor$number_of_variants_turro[all_w_clinvar_no_tumor$Year==x],na.rm=T)}))
clinvar <- rep(unique(all_w_clinvar_no_tumor$Year),times=sapply(unique(all_w_clinvar$Year),function(x){sum(all_w_clinvar_no_tumor$number_of_variants_clinvar[all_w_clinvar_no_tumor$Year==x],na.rm=T)}))
tnamse <- rep(unique(all_w_clinvar_no_tumor$Year),times=sapply(unique(all_w_clinvar$Year),function(x){sum(all_w_clinvar_no_tumor$number_of_variants_tnamse[all_w_clinvar_no_tumor$Year==x],na.rm=T)}))

ks.test(turro,clinvar)
ks.test(turro,tnamse)
  
#tnamse_gene_list<-read_tsv("tnamse_gene_list.txt", col_names = c("gene_list","Freq")) %>% filter(Freq!=0)

#tnamse_gene_list<- tnamse_gene_list %>% arrange(-Freq)
#tnamse_gene_list$Freq<-factor(tnamse_gene_list$Freq, levels=unique(gene_counts$Freq))

#ggplot() + geom_bar(data=tnamse_gene_list, aes(x="gene_list",y=Freq),position="stack", stat="identity", colour="white", fill="black") + 
#geom_bar(data=gene_counts, aes(x="turro",y=Freq),position="stack", stat="identity", colour="white", fill="black") + 
#  theme(legend.position = "none")  + 
#  scale_fill_grey()

