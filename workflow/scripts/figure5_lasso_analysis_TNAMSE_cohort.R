library(tidyverse)
library(ontologyIndex)
library(ontologySimilarity)
library(stats)
library(mlbench)
library(glmnet)
library(plotmo)
library(plotly)
library(pROC)
library(ggrepel)
library(reshape2)
library(directlabels)
library(cowplot)
library(data.table)


## debug
args<-c("results/parse_table/all_cases.tsv",
        "resources/hpo.obo",
        "resources/hpo_categorization_19_12_2022.tsv",
        "results/figure5/",
        "../../")

args = commandArgs(trailingOnly=TRUE)
depth<-args[5]
input_data<-paste0(depth, args[1]) 
input_hpo<-paste0(depth, args[2]) 
input_hpo_categorization<-paste0(depth, args[3]) 
save_path_figures<-paste0(depth, args[4]) 

# original from Hannah
#input_data <- "data_10012022.tsv"
#input_hpo <- "hpo.obo"
#input_hpo_categorization <- "hpo_categorization_08_10_2021.txt"
#save_path_figures <- "figures/"

#set seed

set.seed(1)

#LOAD DATA 

all_cases<-read_tsv(input_data, guess_max = 2000,
                    locale=locale(decimal_mark = ",")) %>% 
  distinct(case_ID_paper,.keep_all = TRUE) %>% 
  filter(!is.na(case_ID_paper))# remove double occuring cases 
#(for security, should already be removed)
#remove double patients (with distinct case_ID_paper)
all_cases_original <- all_cases

# load recent HPO (download 12.02.2021)
hpo<-get_ontology(input_hpo)

#REFORMAT DATA
# format the HPO terms that we can use them as input;remove duplicated HPO terms

all_cases$HPO_Term_IDs <- all_cases$HPO_Term_IDs %>% 
  str_replace_all(";","") %>% strsplit(split=" ") %>%
  sapply(., function(x) unique(x))

#HPO Count per person
all_cases <- all_cases %>%
  mutate(hpo_count=ifelse(is.na(HPO_Term_IDs),0,lengths(HPO_Term_IDs)))

# UNNEST
#first unnest the HPO terms, add a counter as patient ID
expanded_cases <- unnest_longer(all_cases %>% 
                                  mutate(patient_number=1:n()) %>%
                                  group_by(disease_category) %>% 
                                  add_count(name="count_disease_cat") %>% 
                                  ungroup(), 
                                col=HPO_Term_IDs)

#HPO categorization
# Read in Luisas HPO Categories
hpo_categorization<- read_tsv(input_hpo_categorization)
hpo_cat_dict<-hpo_categorization %>% 
  distinct(Subcategory, Category)%>%
  mutate(num=str_split(Category, fixed("-"), simplify = TRUE)[,1])%>%
  mutate(num=as.numeric(num))


hpo_terms_TNAMSE <- unique(hpo_categorization$HPO_term)

# Subcategory is the exact category; Category is the parent category; Reformat:
hpo_categorization01 <- hpo_categorization %>% # reformat exact categories
  group_by(Subcategory) %>% 
  summarise(HPO_terms = list(unique(HPO_term)))

#add descendants of hpo terms
for (sub in 1:nrow(hpo_categorization01)){
  for (term in hpo_categorization01[[2]][[sub]]) {
    descendants<- get_descendants(hpo, term)
    descendants <- descendants[!descendants%in%hpo_terms_TNAMSE]
    hpo_categorization01[[2]][[sub]]=c(hpo_categorization01[[2]][[sub]],
                                       descendants)
  }
  hpo_categorization01[[2]][[sub]]=unique(hpo_categorization01[[2]][[sub]])
}

hpo_categorization02 <- hpo_categorization %>% # reformat parent terms
  group_by(Category) %>% 
  summarise(HPO_terms = list(unique(HPO_term))) %>%
  mutate(Subcategory=paste(Category,"_parental",sep="") ) %>% 
  dplyr::select(-Category)

#add descendants
for (sub in 1:nrow(hpo_categorization02)){
  for (term in hpo_categorization01[[1]][[sub]]) {
    descendants<- get_descendants(hpo, term)
    descendants <- descendants[!descendants%in%hpo_terms_TNAMSE]
    hpo_categorization02[[1]][[sub]]=c(hpo_categorization02[[1]][[sub]],
                                       descendants)
  }
  hpo_categorization02[[1]][[sub]]=unique(hpo_categorization02[[1]][[sub]])
}

hpo_categorization_comb<-rbind(hpo_categorization01,hpo_categorization02) %>% 
  arrange(Subcategory) # combine exact categories and parental categories

# assign patients to HPO categories from categorization
Subcategory<- lapply(hpo_categorization_comb$HPO_terms, 
                   function(x) {
                    return(
                      ifelse(all_cases$case_ID_paper %in%
(expanded_cases[expanded_cases$HPO_Term_IDs %in% unlist(x),]$case_ID_paper),
                                   1, 0)
                          )})

Subcategory<- tibble(as.data.frame(Subcategory), 
                     .name_repair = ~ hpo_categorization_comb$Subcategory) # add column names
all_cases<-cbind(all_cases,Subcategory) # combine with full data set

##Data Cleansing for Lasso

currated_dataset<-all_cases%>%
  mutate(one_for_spread=1)%>% spread(key=disease_category, value=one_for_spread,
                                     fill=0) %>%
  mutate(sequencing_laboratory=as.factor(sequencing_laboratory)) %>% 
  mutate(child=as.integer(adult_child=="child")) %>% 
  mutate(solved=as.integer(solved=="solved")) %>%
  mutate(pedia=as.integer(!is.na(Face2Gene_ID))) %>%
  mutate(sex=as.integer(sex=="female"))


dataset_predict<-currated_dataset %>%
  dplyr::select(case_ID_paper,sex, sequencing_laboratory, child, solved,
                pedia, hpo_count, any_of(hpo_categorization01[[1]])) 

dataset_predict_list<-dataset_predict  %>% as.data.frame()  
# just keep complete cases
dataset_predict_list<-dataset_predict_list[complete.cases(dataset_predict_list),] 

######################### End Data cleaning for lasso ##########################

# Split dataset into train- and testset (validation is done by crossvalidation):
# 80% training, 20% testset
dataset_predict_list[,"train"] <- ifelse(runif(nrow(dataset_predict_list))<0.8,
                                         1,0)
trainset <- dataset_predict_list[dataset_predict_list$train==1,]
testset <- dataset_predict_list[dataset_predict_list$train==0,]
trainColNum <- grep("train",names(trainset))

trainset <- trainset[,-trainColNum]
testset <- testset[,-trainColNum]

### summary of dataset ###
print(paste("In total there are",nrow(all_cases),"cases of which", 
            nrow(dataset_predict_list), 
            "are complete and used for the regression analysis.",
            "The dataset is divided into a training set consisting of",
            nrow(trainset),"cases and a testset consisting of",nrow(testset),
            "cases."))

cases_table=table(dataset_predict_list$train,dataset_predict_list$solved)
rownames(cases_table)=c("test","train")
colnames(cases_table)=c("unsolved","solved")

cases_table

################################################################################
#################### Lasso Regression for HPO subcategories ####################
################## controlling for age, laboratory and pedia ###################
################################################################################

# formula
f_lasso_hpo_subs <-paste0("solved ~ child + sex + sequencing_laboratory + pedia:sequencing_laboratory + '",
                          paste(hpo_categorization01[[1]],collapse="' + '"),"'")

Confounders <-c("child","sex","sequencing_laboratoryB",
                "sequencing_laboratoryC",
                "sequencing_laboratoryD","sequencing_laboratoryE",
                "sequencing_laboratoryA:pedia","sequencing_laboratoryB:pedia")

# define response variable: solved yes/no
y=trainset$solved

#build design matrices
x_lasso_hpo_subs <- as.matrix(cbind(model.matrix(formula("solved ~ child + sex + sequencing_laboratory + pedia:sequencing_laboratory"),
                                                 trainset)[,-1],
                                    trainset[hpo_categorization01[[1]]]))
x_lasso_hpo_subs <- x_lasso_hpo_subs[,!colnames(x_lasso_hpo_subs) %in% c("sequencing_laboratoryC:pedia","sequencing_laboratoryD:pedia","sequencing_laboratoryE:pedia")]
x_lasso_hpo_subs_test <- as.matrix(cbind(model.matrix(formula("solved ~ child + sex + sequencing_laboratory + pedia:sequencing_laboratory"),
                                                      testset)[,-1],
                                         testset[hpo_categorization01[[1]]]))
x_lasso_hpo_subs_test <- x_lasso_hpo_subs_test[,!colnames(x_lasso_hpo_subs_test) %in% c("sequencing_laboratoryC:pedia","sequencing_laboratoryD:pedia","sequencing_laboratoryE:pedia")]


# confounder should not be penalized in the model
pen_fact<-as.integer(!colnames(x_lasso_hpo_subs) %in% Confounders)

# fit lasso with maximizing auc via 10-fold crossvalidation
lasso_hpo_subs <- cv.glmnet(x_lasso_hpo_subs, y, alpha=1, family="binomial",
                            penalty.factor=pen_fact, type.measure = "auc",
                            nfolds=10)

# plot coefficient paths seperately

p=list()
for(i in 1:12){
  relevant_cats<-hpo_cat_dict %>% 
    filter(num==i)
  
  Koeff=as.data.frame(cbind(lasso_hpo_subs$lambda,t(as.matrix(lasso_hpo_subs$glmnet.fit$beta))))
  Koeff=Koeff[,colnames(Koeff)=="V1"| colnames(Koeff) %in% relevant_cats$Subcategory]
   df_long=reshape2::melt(Koeff,id.vars="V1")
  p[[i]] <- ggplot(df_long,aes(x=log(V1),y=value,color=variable))+geom_line() +
    labs(x=expression(ln(lambda)),y="Coefficients") +
    geom_vline(xintercept=log(lasso_hpo_subs$lambda.min), color= "black") +
    scale_colour_discrete(guide = 'none') +
    scale_x_reverse(expand=expansion(add=c(0,10))) +
    geom_text_repel(data = . %>% filter(V1==min(V1)),
                    mapping=aes(label = paste0("  ",variable),segment.color="grey",hjust=0),
                    nudge_x =0.5,
                    direction = "y") +
    theme_minimal() 
  
    ggsave(filename=paste0(save_path_figures,"Koeffpfad_Lasso_",i,".pdf"),plot=p[[i]], width=7,height=5)
}

print(plot_grid(plotlist=p,ncol=2))

table_hpo_subs=data.table(subgroup=hpo_categorization01[[1]],count=rep(0,length(hpo_categorization01[[1]])))
Koeff=as.data.frame(cbind(lasso_hpo_subs$lambda,t(as.matrix(lasso_hpo_subs$glmnet.fit$beta))))

for(i in 1:nrow(table_hpo_subs)){
  table_hpo_subs$count[i] = sum(trainset[hpo_categorization01[[1]][i]])
  table_hpo_subs$coeff_at_best[i]=Koeff[Koeff$V1==lasso_hpo_subs$lambda.min,table_hpo_subs$subgroup[i]]
}

## choose coefficients that are influential and relevant for at least 5% of the trainset
table_hpo_subs$include = table_hpo_subs$coeff_at_best!=0 & table_hpo_subs$count>=0.05*nrow(trainset)
#table_hpo_subs

Koeff=Koeff[,colnames(Koeff)=="V1"|colnames(Koeff)%in% table_hpo_subs$subgroup[table_hpo_subs$include]]
colnames(Koeff)=str_replace(colnames(Koeff),"Abnormality of","Ao")
colnames(Koeff)=str_replace(colnames(Koeff),"Only ","")
df_long=reshape2::melt(Koeff,id.vars="V1")
df_long<-df_long %>% left_join(hpo_cat_dict, by=c("variable"="Subcategory"))

plot_manusscript <- ggplot(df_long %>% filter(log(V1)>-9),
                           aes(x=log(V1),y=value,group=variable, color=Category))+
  geom_line() +
  labs(x=expression(ln(lambda)),y="Coefficients") +
  geom_vline(xintercept=log(lasso_hpo_subs$lambda.min), color= "black") +
  scale_colour_discrete(guide = 'none') +
  scale_x_reverse(expand=expansion(add=c(0,6))) +
  geom_text_repel(data = . %>% filter(V1==min(V1)),
                  mapping=aes(label = paste0("  ",variable),segment.color="grey",
                              hjust=0),
                  nudge_x =0.5,
                  direction = "y") +
  theme_minimal()

  ggsave(filename=paste0(save_path_figures,"lasso_coefficients_manusscript.pdf"),plot=plot_manusscript,width=7,
         height=5)

plot_manusscript

print(paste0("Lambda chosen by 10-fold crossvalidation: ",
             expression(lambda),"=",lasso_hpo_subs$lambda.min))
Ausgabe=cbind(coef(lasso_hpo_subs,s=lasso_hpo_subs$lambda.min),
              exp(coef(lasso_hpo_subs,s=lasso_hpo_subs$lambda.min)))
colnames(Ausgabe)=c("Coefficients","Odds Ratios")
print(c("The corresponding coefficients and odds ratios are given by: ", 
        Ausgabe))
print(paste0("Number of non-zero coefficients: ",
             sum(coef(lasso_hpo_subs, s=lasso_hpo_subs$lambda.min)!=0)))

lasso_hpo_subs_fitted= predict.glmnet(lasso_hpo_subs$glmnet.fit,
                                      x_lasso_hpo_subs_test,
                                      s=lasso_hpo_subs$lambda.min)
roc(testset$solved, as.numeric(lasso_hpo_subs_fitted) ,
    percent=F, stratified=FALSE, plot=TRUE, 
    grid=TRUE, show.thres=TRUE, legacy.axes = TRUE, 
    reuse.auc = TRUE,print.auc = TRUE, print.thres.col = "blue", ci=TRUE,
    ci.type="bars", print.thres.cex = 0.7, 
    main = paste("ROC curve using","(N = ",length(lasso_hpo_subs_fitted),")"))

testset_outcome<-tibble(
  solved=testset$solved,
  predicted=1/(1+exp(-lasso_hpo_subs_fitted))
)

ggplot(testset_outcome, aes(x=predicted, color=solved==1))+
  geom_histogram(aes(fill=solved==1), position="identity", alpha=0.5)

ggplot(testset_outcome, aes(x=predicted, fill=solved==1))+
  geom_density(alpha=0.5)


########### test using the mean value of confounders in training set ###########
# set confounder columns to mean of trainset

x_lasso_hpo_subs_test_no_confounder <- x_lasso_hpo_subs_test
for(c in Confounders){
  x_lasso_hpo_subs_test_no_confounder[,c] = mean(x_lasso_hpo_subs[,c])
}

lasso_hpo_subs_fitted_no_confounders=
  predict.glmnet(lasso_hpo_subs$glmnet.fit, x_lasso_hpo_subs_test_no_confounder,
                 s=lasso_hpo_subs$lambda.min)
roc(testset$solved, as.numeric(lasso_hpo_subs_fitted_no_confounders) , percent=F, 
    stratified=FALSE, plot=TRUE, grid=TRUE, show.thres=TRUE, legacy.axes = TRUE, 
    reuse.auc = TRUE,print.auc = TRUE, print.thres.col = "blue", ci=TRUE, 
    ci.type="bars", print.thres.cex = 0.7, 
    main = paste("ROC curve using","(N = ",length(lasso_hpo_subs_fitted),")"))





