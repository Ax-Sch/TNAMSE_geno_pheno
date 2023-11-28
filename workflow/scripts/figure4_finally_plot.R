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

# Check how many variants are in tumor/ACMGv2 genes:
print("Number of variants in tumor genes in TNAMSE, Turro and ClinVar + histogram for ClinVar:")
sum((all_w_clinvar %>% filter(cancer_gene==TRUE))$number_of_variants_tnamse, na.rm=TRUE)
sum((all_w_clinvar %>% filter(cancer_gene==TRUE))$number_of_variants_turro, na.rm=TRUE)
cv_tumor<-(all_w_clinvar %>% filter(cancer_gene==TRUE))$number_of_variants_clinvar
sum(cv_tumor, na.rm=TRUE)
qplot(cv_tumor)+ theme_minimal()

# Further prepare the data
all_w_clinvar_no_tumor<-all_w_clinvar %>%
  mutate(number_of_variants_clinvar=replace_na(data=number_of_variants_clinvar, replace=0))%>%
  filter(cancer_gene==FALSE)  %>%
  arrange(-number_of_variants_clinvar)%>% 
  mutate(cumsum_clinvar=cumsum(number_of_variants_clinvar))

total_clinvar<-sum(all_w_clinvar_no_tumor$number_of_variants_clinvar, na.rm=TRUE)

max_TNAMSE<-max(all_w_clinvar_no_tumor$number_of_variants_tnamse, na.rm = TRUE)
max_Turro<-max(all_w_clinvar_no_tumor$number_of_variants_turro, na.rm = TRUE)
max_clinvar<-max(all_w_clinvar_no_tumor$number_of_variants_clinvar, na.rm=TRUE)


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
                        labels=c("-2000","2001-2010","2011-2022")))

# change years / clinvar quartile to factor
all_w_clinvar_no_tumor<-all_w_clinvar_no_tumor %>% 
  arrange(str_sort(year_range, numeric=TRUE))%>%
  mutate(year_range=factor(year_range, levels=unique(year_range)))%>% 
  arrange((clinvar_quarter))%>%
  mutate(clinvar_quarter=factor(clinvar_quarter, levels=unique(clinvar_quarter)))

write_tsv(x=all_w_clinvar_no_tumor, 
          file="all_w_clinvar_no_tumor.tsv")


print("genes in quarters")
print(table(all_w_clinvar_no_tumor$clinvar_quarter))

cv_quarters<-all_w_clinvar_no_tumor%>% 
  group_by(clinvar_quarter)%>%
  summarise(max_cv=max(cumsum_clinvar))%>%
  distinct(clinvar_quarter, max_cv)

print("cumulative vars relative to quarters")
print(cv_quarters)


# Scatter plots with ClinVar as reference
freq_clinvar_vs_tnamse<-ggplot(all_w_clinvar_no_tumor)+
  geom_line(aes(x=cumsum_clinvar,y=number_of_variants_clinvar/max_clinvar), size=0.5)+
  geom_point(aes(x=cumsum_clinvar,y=number_of_variants_clinvar/max_clinvar), size=1)+
  geom_jitter(data=all_w_clinvar_no_tumor[!is.na(all_w_clinvar_no_tumor$number_of_variants_tnamse),], size=1.0, alpha=0.7, width=0.05,height=0.05, 
              aes(x=cumsum_clinvar,y=number_of_variants_tnamse/max_TNAMSE, color=year_range))+
  theme_bw()
print(freq_clinvar_vs_tnamse)
ggsave(freq_clinvar_vs_tnamse, filename = "freq_clinvar_vs_tnamse.pdf", width=4.5, height=3.0)

freq_clinvar_vs_turro<-ggplot(all_w_clinvar_no_tumor)+
  geom_line(aes(x=cumsum_clinvar,y=number_of_variants_clinvar/max_clinvar), size=0.5)+
  geom_point(aes(x=cumsum_clinvar,y=number_of_variants_clinvar/max_clinvar), size=1)+
  geom_jitter(data=all_w_clinvar_no_tumor[!is.na(all_w_clinvar_no_tumor$number_of_variants_turro),], size=1.0, alpha=0.7, width=0.01,height=0.01, 
              aes(x=cumsum_clinvar,y=number_of_variants_turro/max_Turro, color=year_range))+
  theme_bw()
print(freq_clinvar_vs_turro)
ggsave(freq_clinvar_vs_turro, filename = "freq_clinvar_vs_turro.pdf", width=4.5, height=3.0)



# Bar plot number of gene occurance 
plot_number<-ggplot() + 
  geom_bar(data=all_w_clinvar_no_tumor %>% arrange(-number_of_variants_turro), aes(x="Turro et al.",y=number_of_variants_turro), color="grey40",position="stack", fill="white",stat="identity", width=0.7) + 
  geom_bar(data=all_w_clinvar_no_tumor %>% arrange(-number_of_variants_tnamse), aes(x="TNAMSE",y=number_of_variants_tnamse), color="grey40",position="stack", fill="white",stat="identity", width=0.7) + 
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45))
print(plot_number)
ggsave(file="plot_numbers.pdf",plot_number, width=1.8, height=3.2)


# Years first report
plot_year_tn<-ggplot(data=all_w_clinvar_no_tumor %>% arrange(-number_of_variants_tnamse)) + # %>% filter(!is.na(year_range))
  geom_bar( aes(x=year_range,y=number_of_variants_tnamse, color=year_range),position="stack", fill="white",stat="identity", width=0.7) + 
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, hjust=0.95))
print(plot_year_tn)
ggsave(file="plot_year_tn.pdf",plot_year_tn, width=4.0, height=3.2)

plot_year_turro<-ggplot(data=all_w_clinvar_no_tumor %>% arrange(-number_of_variants_turro) ) + #%>% filter(!is.na(year_range))
  geom_bar( aes(x=year_range,y=number_of_variants_turro, color=year_range),position="stack", fill="white",stat="identity", width=0.7) + 
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, hjust=0.95))
print(plot_year_turro)
ggsave(file="plot_year_turro.pdf",plot_year_turro, width=4.0, height=3.2)


# ClinVar Quartiles
plot_clinvar_clinvar<-ggplot() + 
  geom_bar(data=all_w_clinvar_no_tumor %>% arrange(-number_of_variants_clinvar), 
           aes(x=clinvar_quarter,y=number_of_variants_clinvar, color=clinvar_quarter),position="stack", fill="white",stat="identity") + 
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45))
print(plot_clinvar_clinvar)
ggsave(file="plot_clinvar_clinvar.pdf",plot_clinvar_clinvar, width=4.5, height=7)

plot_clinvar_turro<-ggplot() + 
  geom_bar(data=all_w_clinvar_no_tumor %>% arrange(-number_of_variants_turro), 
           aes(x=clinvar_quarter,y=number_of_variants_turro, color=clinvar_quarter),position="stack", fill="white",stat="identity", width=0.7) + 
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, hjust=0.95))
print(plot_clinvar_turro)
ggsave(file="plot_clinvar_turro.pdf",plot_clinvar_turro, width=4.0, height=3.2)

plot_clinvar_tn<-ggplot() + 
  geom_bar(data=all_w_clinvar_no_tumor %>% arrange(-number_of_variants_tnamse), 
           aes(x=clinvar_quarter,y=number_of_variants_tnamse, color=clinvar_quarter),position="stack", fill="white",stat="identity", width=0.7) + 
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, hjust=0.95))
print(plot_clinvar_tn)
ggsave(file="plot_clinvar_tn.pdf",plot_clinvar_tn, width=4.0, height=3.2)


#### Kolmogorov-Smirnoff-Test
turro <- rep(unique(all_w_clinvar_no_tumor$Year),times=
            sapply(unique(all_w_clinvar$Year),
                   function(x){sum(all_w_clinvar_no_tumor$number_of_variants_turro[all_w_clinvar_no_tumor$Year==x],na.rm=T)}
                   ))
clinvar <- rep(unique(all_w_clinvar_no_tumor$Year),times=
              sapply(unique(all_w_clinvar$Year),
                   function(x){sum(all_w_clinvar_no_tumor$number_of_variants_clinvar[all_w_clinvar_no_tumor$Year==x],na.rm=T)}
                   ))
tnamse <- rep(unique(all_w_clinvar_no_tumor$Year),times=
              sapply(unique(all_w_clinvar$Year),
                   function(x){sum(all_w_clinvar_no_tumor$number_of_variants_tnamse[all_w_clinvar_no_tumor$Year==x],na.rm=T)}
                   ))

ks.test(turro,clinvar)
ks.test(turro,tnamse)

