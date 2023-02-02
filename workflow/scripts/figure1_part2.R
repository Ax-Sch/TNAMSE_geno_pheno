### Master ontology created
library(ontologyIndex)
library(Rcpp)
library(ontologySimilarity)
library(tidyverse)
library(umap)
library(ggrepel)
library(flexclust)

# for debugging; if the last parameter is equal to "redo", the similarity matrix (time consuming!, ca. 6h) will be recalculated
args=c("results/parse_table/all_cases_wHighEvNovel.tsv",
       "resources/hpo.obo", 
       "resources/genes_to_phenotype.txt",
       "results/figure4/Turro_variant_table.txt",
       "../../",
       "redo")

args = commandArgs(trailingOnly=TRUE)
depth<-args[5]
in_file<-paste0(depth, args[1]) #
hpo_obo<-paste0(depth, args[2]) #"ressources/hpo.obo"
gene_to_pheno_path<-paste0(depth, args[3]) #"ressources/hpo_categorization_luisa_01_07_2021.txt"
turro_path<-paste0(depth, args[4])
redo=args[6]

hpo<-get_ontology(hpo_obo)
blacklist_hpos = c("HP:0000006", "HP:0000007", "HP:0001417", "HP:0001419", "HP:0001423", "HP:0001428", "HP:0001450", "HP:0040284", "HP:0040283")
only_phenotypes_hpo = get_descendants(hpo, "HP:0000118", exclude_roots = FALSE)
hpos_to_keep =only_phenotypes_hpo[!only_phenotypes_hpo %in% blacklist_hpos]


#LOAD DATA 
TNAMSE_data<-read_tsv(in_file, guess_max = 2000, locale=locale(decimal_mark = ",")) %>% 
  distinct(case_ID_paper,.keep_all = TRUE) # remove double occuring cases (for security, should already be removed)
# format the HPO terms that we can use them as input; remove duplicated HPO terms
gene_to_pheno <- read.table(gene_to_pheno_path, sep="\t", quote="", stringsAsFactors=FALSE, header=TRUE)
turro<-read_tsv(turro_path)

TNAMSE_data <- TNAMSE_data %>% filter(is.na(HPO_Term_IDs)==FALSE)

TNAMSE_data$HPO_Term_IDs <- TNAMSE_data$HPO_Term_IDs %>% 
  str_replace_all(";","") %>% strsplit(split=" ") %>%
  sapply(., function(x) unique(x)) 

TNAMSE_data<-TNAMSE_data %>% unnest_longer(HPO_Term_IDs)
TNAMSE_data <- TNAMSE_data[TNAMSE_data$HPO_Term_IDs %in% hpos_to_keep,] # remove HPO terms not wanted

TNAMSE_data_red <- TNAMSE_data %>% distinct(case_ID_paper, .keep_all = TRUE) %>% select(-HPO_Term_IDs)

"HPO_Term_IDs"
for (x in 1:nrow(TNAMSE_data_red)) {
  TNAMSE_data_red$HPO_term_IDs[[x]] <- (TNAMSE_data[which(TNAMSE_data$case_ID_paper == TNAMSE_data_red[x,]$case_ID_paper), ]$HPO_Term_IDs)
}

# downloaded 10.04.2021; manually modified column header
gene_to_pheno = gene_to_pheno[gene_to_pheno$HPO_Term_ID %in% hpos_to_keep, ]
overall_disease = unique(gene_to_pheno$disease_ID_for_link)


list_of_phenotypes_HPO<-data.frame()
for (disease in overall_disease) {
  gene_to_pheno_disease<-gene_to_pheno[gene_to_pheno$disease_ID_for_link == disease,]
  list_of_phenotypes_HPO<- rbind(list_of_phenotypes_HPO,
                                 data.frame(
                                   case_ID_paper=disease, 
                                   HPO_term_IDs=I(list((gene_to_pheno_disease$HPO_Term_ID))),
                                   disease_category="HPO",
                                   sequencing_laboratory=substring(disease,1,4),
                                   Disease_gene=I(list(unique(gene_to_pheno_disease$entrez_gene_symbol)))
                                 )
  )
}

print(paste("number of patients in similarity matrix:", nrow(TNAMSE_data_red)))
print(paste("number of diseases in similarity matrix:", nrow(list_of_phenotypes_HPO)))
print(paste("number of diseases in similarity matrix, overall_diseases:", nrow(overall_disease), length(overall_disease)))

library(plyr)
TNAMSE_and_HPO<-rbind.fill(
  TNAMSE_data_red,
  list_of_phenotypes_HPO
)
detach("package:plyr", unload=TRUE)

information_content <- descendants_IC(hpo)

if (redo == "redo"){
  print("repeating similarity grid")
  master_sim_mat <- get_sim_grid(ontology=hpo, information_content=information_content,  term_sets=TNAMSE_and_HPO$HPO_term_IDs, term_sim_method = "resnik", combine="average")
  write_rds(x = master_sim_mat, file = "master_sim_mat.RDS")
}
#master_sim_mat<-read_rds(file="../../results/figure1/master_sim_mat.RDS")
master_sim_mat<-read_rds(file="master_sim_mat.RDS")

set.seed(1)
custom.settings = umap.defaults
custom.settings$input = "dist"
custom.settings$n_components = 4

#why **2 ?
res_umap <- umap(as.matrix(max(master_sim_mat) - master_sim_mat) ** 2, config=custom.settings)

colnames(res_umap$layout)<-paste0("dim", 1:(custom.settings$n_components))
TNAMSE_and_HPO <- cbind(TNAMSE_and_HPO,res_umap$layout )

TNAMSE_and_HPO<- TNAMSE_and_HPO%>%
  mutate(dim1=-dim1,
         dim2=-dim2)

ggplot(TNAMSE_and_HPO, aes(x=dim1, y=dim2, color=disease_category)) +
  geom_point(alpha=0.2)


number_of_clusters = 80
set.seed(1)
only_HPO<-TNAMSE_and_HPO %>% filter(disease_category == "HPO")

clusters = kcca(only_HPO %>% dplyr::select(one_of(paste0("dim", 1:(custom.settings$n_components)))), k=number_of_clusters, kccaFamily("kmeans"))
#plot(result_only_disease[,1], result_only_disease[,2], col=rainbow(max(clusters$cluster))[clusters$cluster+1], pch=20)



only_HPO$cluster<-predict(clusters)
TNAMSE_and_HPO$cluster_pred<-predict(clusters, newdata=TNAMSE_and_HPO %>% dplyr::select(one_of(paste0("dim", 1:(custom.settings$n_components)))))

most_frequent_hpos_per_cluster<-only_HPO %>% group_by(cluster) %>% 
  add_count(name="count_of_patients_in_cluster") %>%
  ungroup() %>%
  unnest_longer(col=HPO_term_IDs) %>% 
  distinct(HPO_term_IDs,case_ID_paper, .keep_all = TRUE) %>%
  group_by(cluster,HPO_term_IDs) %>%
  add_count(name="count_of_HPO_term_in_Cluster")%>%
  distinct(cluster,HPO_term_IDs,count_of_HPO_term_in_Cluster,count_of_patients_in_cluster) %>%
  arrange(cluster,-count_of_HPO_term_in_Cluster)%>%
  mutate(proportion_w_HPO=round(count_of_HPO_term_in_Cluster/count_of_patients_in_cluster, 2))%>%
  group_by(cluster)%>%
  slice_max(order_by = count_of_HPO_term_in_Cluster, n = 5, with_ties=FALSE)

most_frequent_hpos_per_cluster$HPO_Description<-sapply(most_frequent_hpos_per_cluster$HPO_term_IDs, function(x) hpo$name[hpo$id==x])

cluster_descriptions<-most_frequent_hpos_per_cluster %>% 
  group_by(cluster, count_of_patients_in_cluster) %>%
  mutate(description=paste0(paste(proportion_w_HPO,HPO_Description), collapse = "\n")) %>%
  distinct(cluster,count_of_patients_in_cluster,description)

cluster_stats_TNAMSE_cohort<-TNAMSE_and_HPO %>% 
  dplyr::filter(disease_category!="HPO" | is.na(disease_category)) %>%
  group_by(cluster_pred)%>%
  add_count(name="cluster_size")%>%
  group_by(cluster_pred,cluster_size)%>%
  mutate(is_solved=mean(solved=="solved", na.rm=TRUE)) %>%
  distinct(is_solved,cluster_pred,cluster_size)

cluster_descriptions<- cbind(as.data.frame(cluster_descriptions), clusters@centers)
write_tsv(x=cluster_descriptions, "cluster_descriptions.tsv")

#clusters_for_diseases <- read.table("clusters_with_clinical_annotation.txt", sep="\t", header=T)
#cluster_descriptions$manually_annotated_category<-clusters_for_diseases$Category


TNAMSE_and_HPO <-TNAMSE_and_HPO %>% 
  left_join(only_HPO %>% dplyr::select(case_ID_paper, cluster),
            by=c("case_ID_paper"="case_ID_paper"), )%>%
  left_join(cluster_stats_TNAMSE_cohort,
            by=c("cluster_pred"="cluster_pred"), )






plot_interim<-ggplot()+ theme_minimal() + theme(legend.position="bottom",panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_point(data = TNAMSE_and_HPO %>% filter(disease_category=="HPO"),
             aes(x = dim1, y = dim2), alpha=0.5, color="lightgrey") +
  geom_point(data = TNAMSE_and_HPO %>% filter(disease_category!="HPO" & !novel_disease_gene),
             aes(x = dim1, y = dim2, color = disease_category), size=3) + 
  geom_point(data = TNAMSE_and_HPO %>% filter(disease_category!="HPO" & novel_disease_gene),
             aes(x = dim1, y = dim2, color = disease_category), shape=17, size=5) + 
  geom_label_repel(data = cluster_descriptions, aes(x=dim1, y=dim2, label=description), size=2, color="black", show.legend = FALSE)

ggsave("plot_w_label.pdf", plot = plot_interim, width=50, height=50, units="cm", dpi=500, useDingbats=FALSE)



plot_interim<-ggplot()+ theme_minimal() + theme(legend.position="bottom",panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_point(data = TNAMSE_and_HPO %>% filter(disease_category=="HPO"),
             aes(x = dim1, y = dim2), alpha=0.5, color="lightgrey") +
  geom_point(data = TNAMSE_and_HPO %>% filter(disease_category!="HPO" & !novel_disease_gene),
             aes(x = dim1, y = dim2, color = disease_category), size=3, shape=19) + 
  geom_point(data = TNAMSE_and_HPO %>% filter(disease_category!="HPO" & novel_disease_gene),
             aes(x = dim1, y = dim2, color = disease_category), fill="black", shape=24, size=4.0, stroke=1.5)

plot_interim

ggsave("plot_wo_label.pdf", plot = plot_interim, width=25, height=25, units="cm", dpi=500, useDingbats=FALSE)


# code taken from:
#https://stackoverflow.com/questions/28521145/r-calculate-and-plot-difference-between-two-density-countours
# 12.04.2021

library(MASS)
library(reshape2) # For melt function

density_plot<-function(dataset1, dataset2,title_par="no title given"){
  # Calculate the common x and y range for geyser1 and geyser2
  xrng = range(c(dataset1$dim1, dataset2$dim1))
  xrng = range(c(xrng+1, xrng-1))
  yrng =  range(c(dataset1$dim2, dataset2$dim2))
  yrng = range(c(yrng+1, yrng-1))
  
  # Calculate the 2d density estimate over the common range
  d1 = kde2d(dataset1$dim1,dataset1$dim2, lims=c(xrng, yrng), n=200)
  d2 = kde2d(dataset2$dim1, dataset2$dim2, lims=c(xrng, yrng), n=200)
  
  # Confirm that the grid points for each density estimate are identical
  identical(d1$x, d2$x) # TRUE
  identical(d1$y, d2$y) # TRUE
  
  # Calculate the difference between the 2d density estimates
  diff12 = d1 
  diff12$z = d2$z - d1$z
  
  ## Melt data into long format
  # First, add row and column names (x and y grid values) to the z-value matrix
  rownames(diff12$z) = diff12$x
  colnames(diff12$z) = diff12$y
  
  # Now melt it to long format
  diff12.m = melt(diff12$z, id.var=rownames(diff12))
  names(diff12.m) = c("dim1","dim2","z")
  
  # Plot difference between geyser2 and geyser1 density
  return(ggplot(diff12.m, aes(dim1, dim2, z=z, fill=z)) +
           ggtitle(title_par) +
           stat_contour(aes(colour=..level..), binwidth=0.001) +
           scale_fill_gradient2(low="red",mid="white", high="blue", midpoint=0) +
           scale_colour_gradient2(low="red", mid="white", high="blue", midpoint=0) +
           geom_point(data = TNAMSE_and_HPO %>% filter(disease_category=="HPO"),
                      aes(x = dim1, y = dim2), alpha=0.15, color="lightgrey",inherit.aes = FALSE)+
           geom_point(data = TNAMSE_and_HPO %>% filter(disease_category!="HPO"),
                      aes(x = dim1, y = dim2), size=3, inherit.aes = FALSE, alpha=0.07,color="green")  + theme_bw())
  
}

solved_TN<-TNAMSE_and_HPO %>% filter(disease_category!="HPO" | is.na(disease_category) ) %>% filter(solved=="solved")
unsolved_TN<-TNAMSE_and_HPO %>% filter(disease_category!="HPO" | is.na(disease_category)) %>% filter(solved!="solved")

ggsave(density_plot(solved_TN, unsolved_TN, "solved vs unsolved"), file="solved_unsolved.pdf", useDingbats=FALSE)

TN<-TNAMSE_and_HPO %>% filter(disease_category!="HPO" | is.na(disease_category))
HPO<-TNAMSE_and_HPO %>% filter(disease_category=="HPO")

ggsave(density_plot(TN, HPO, "distribution HPO vs TNAMSE"), file="HPO_vs_NAMSE.pdf", useDingbats=FALSE)


lb_TB<-TN %>% filter(sequencing_laboratory=="E")
lb_other<-TN %>% filter(sequencing_laboratory!="E")

ggsave(density_plot(lb_TB, lb_other, "sequenced in E"), file="sequenced_in_E.pdf", useDingbats=FALSE)


lb_BN<-TN %>% filter(sequencing_laboratory=="B")
lb_other<-TN %>% filter(sequencing_laboratory!="B")

ggsave(density_plot(lb_BN, lb_other, "sequenced in B"), file="sequenced_in_B.pdf", useDingbats=FALSE)


lb_MN<-TN %>% filter(sequencing_laboratory=="D")
lb_other<-TN %>% filter(sequencing_laboratory!="D")

ggsave(density_plot(lb_MN, lb_other, "sequenced in D"), file="sequenced_in_D.pdf", useDingbats=FALSE)


lb_MN<-TN %>%filter(sequencing_laboratory=="C")
lb_other<-TN %>% filter(sequencing_laboratory!="C")

ggsave(density_plot(lb_MN, lb_other, "sequenced in LMU"), file="sequenced_in_LMU.pdf", useDingbats=FALSE)


lb_kind<-TN %>% filter(adult_child=="child")
lb_erw<-TN %>%  filter(adult_child!="child")

ggsave(density_plot(lb_kind, lb_erw, "Child vs. Adult"), file="adult_child.pdf", useDingbats=FALSE)
