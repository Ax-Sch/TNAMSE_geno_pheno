library(data.table)
library(tidyverse)
library(readxl)
library(cowplot)
library(ggpubr)
library(boot)

data <- fread("data_pedia.tsv")
save_path_figures <- "figures/"

###################### analysis of top pedia score #############################

################ Plots + Univariate logistic regression#########################

#create plots
comparisons = list(c("1 - solved and modelled","2 - solved but not modelled"),
                   c("1 - solved and modelled","3 - unsolved"),
                   c("1 - solved and modelled","4 - unclear"))

vs_unsolved <- ggplot(data %>% filter(group %in% c("1 - solved and modelled","3 - unsolved")),aes(x=pedia_top,color=group))+
  geom_density(aes(fill=group),alpha=0.5)
vs_unclear <- ggplot(data %>% filter(group %in% c("1 - solved and modelled","4 - unclear")),aes(x=pedia_top,color=group))+
  geom_density(aes(fill=group),alpha=0.5)
vs_not_modelled <- ggplot(data %>% filter(group %in% c("1 - solved and modelled","2 - solved but not modelled")),aes(x=pedia_top,color=group))+
  geom_density(aes(fill=group),alpha=0.5)

plot_grid(plotlist=list(vs_not_modelled,vs_unsolved,vs_unclear),ncol=2)

boxplot_plot <- ggplot(data,aes(x=group,y=pedia_top,color=group))+
  geom_boxplot(aes(fill=group),alpha=0.5) +
  stat_compare_means(method = "t.test", label="p.signif", comparisons = comparisons,
                     method.args = list(alternative = "greater", var.equal=F))
boxplot_plot

### linear regression: top pedia score ~ group (categorial effect, reference: solved and modelled)
summary(glm(pedia_top ~ group, data= data, family="gaussian"))
##### result: the pedia score differs significantly between each group and the reference group 1,
##### The estimated effect is negative, i.e. group 1 has a significantly higher top pedia score

### logistic regression: solved and modelled/not(solved and modelled) ~ top pedia score
summary(glm(group=="1 - solved and modelled"~pedia_top,data ,family="binomial"))
##### result: A higher pedia score is significantly associated with a higher prob. of being in group 1


############### Analysis of accuracies: pedia vs. other approaches #############

############ functions to calculate accuracies ################
accuracies_stat_pedia<- function(data,indices){
  d=data[indices,]
  
  tops=numeric(100)
  for(l in 1:100){
    tops[l]=sum(d[[paste0("top",l,"_pedia")]])/nrow(d)
  }
  
  return(tops)
}

accuracies_stat_cadd<- function(data,indices){
  d=data[indices,]
  
  tops=numeric(100)
  for(l in 1:100){
    tops[l]=sum(d[[paste0("top",l,"_cadd")]])/nrow(d)
  }
  
  return(tops)
}

accuracies_stat_gestalt<- function(data,indices){
  d=data[indices,]
  
  tops=numeric(100)
  for(l in 1:100){
    tops[l]=sum(d[[paste0("top",l,"_gestalt")]])/nrow(d)
  }
  
  return(tops)
}

accuracies_stat_phenomizer<- function(data,indices){
  d=data[indices,]
  
  tops=numeric(100)
  for(l in 1:100){
    tops[l]=sum(d[[paste0("top",l,"_phenomizer")]])/nrow(d)
  }
  
  return(tops)
}

accuracies_stat_boqa<- function(data,indices){
  d=data[indices,]
  
  tops=numeric(100)
  for(l in 1:100){
    tops[l]=sum(d[[paste0("top",l,"_boqa")]])/nrow(d)
  }
  
  return(tops)
}

##### filter for group 1 and bootstrapping
data_accuracies <- data %>% filter(group=="1 - solved and modelled")

set.seed(1)
boot_pedia <- boot(data_accuracies,accuracies_stat_pedia,10000)

set.seed(1)
boot_cadd <- boot(data_accuracies,accuracies_stat_cadd,10000)

set.seed(1)
boot_gestalt <- boot(data_accuracies,accuracies_stat_gestalt,10000)

set.seed(1)
boot_phenomizer <- boot(data_accuracies,accuracies_stat_phenomizer,10000)

set.seed(1)
boot_boqa <- boot(data_accuracies,accuracies_stat_boqa,10000)

mean_PEDIA = boot_pedia$t0
mean_CADD = boot_cadd$t0
mean_DeepGestalt = boot_gestalt$t0
mean_Phenomizer = boot_phenomizer$t0
mean_Boqa = boot_boqa$t0

lower_PEDIA <- lower_CADD <- lower_DeepGestalt <- lower_Phenomizer<- lower_Boqa <-numeric(100)
upper_PEDIA <- upper_CADD <- upper_DeepGestalt <- upper_Phenomizer<- upper_Boqa <- numeric(100)

for(l in 1:100){
  lower_PEDIA[l] = if(mean_PEDIA[l]!=1){boot.ci(boot_pedia,type="bca",index=l)$bca[4]}else{1}
  lower_CADD[l] = if(mean_CADD[l]!=1){boot.ci(boot_cadd,type="bca",index=l)$bca[4]}else{1}
  lower_DeepGestalt[l] = if(mean_DeepGestalt[l]!=1){boot.ci(boot_gestalt,type="bca",index=l)$bca[4]}else{1}
  lower_Phenomizer[l] = if(mean_Phenomizer[l]!=1){boot.ci(boot_phenomizer,type="bca",index=l)$bca[4]}else{1}
  lower_Boqa[l] = if(mean_Boqa[l]!=1&mean_Boqa[l]!=0){boot.ci(boot_boqa,type="bca",index=l)$bca[4]}else{mean_Boqa[l]}
  
  upper_PEDIA[l] = if(mean_PEDIA[l]!=1){boot.ci(boot_pedia,type="bca",index=l)$bca[5]}else{1}
  upper_CADD[l] = if(mean_CADD[l]!=1){boot.ci(boot_cadd,type="bca",index=l)$bca[5]}else{1}
  upper_DeepGestalt[l] = if(mean_DeepGestalt[l]!=1){boot.ci(boot_gestalt,type="bca",index=l)$bca[5]}else{1}
  upper_Phenomizer[l] = if(mean_Phenomizer[l]!=1){boot.ci(boot_phenomizer,type="bca",index=l)$bca[5]}else{1}
  upper_Boqa[l] = if(mean_Boqa[l]!=1&mean_Boqa[l]!=0){boot.ci(boot_boqa,type="bca",index=l)$bca[5]}else{mean_Boqa[l]}
}


accuracies <- data.table(x=rep(1:100,times=3),
                         type=rep(c("mean","lower","upper"),each=100),
                         PEDIA=c(mean_PEDIA,lower_PEDIA,upper_PEDIA),
                         CADD=c(mean_CADD,lower_CADD,upper_CADD),
                         DeepGestalt=c(mean_DeepGestalt,lower_DeepGestalt,upper_DeepGestalt),
                         Phenomizer=c(mean_Phenomizer,lower_Phenomizer,upper_Phenomizer),
                         Boqa=c(mean_Boqa,lower_Boqa,upper_Boqa)) %>% 
  pivot_longer(cols=c("PEDIA","CADD","DeepGestalt","Phenomizer","Boqa")) %>% rename("Scoring approach"=name) %>%
  pivot_wider(names_from = 2)

ggplot(accuracies,aes(x=x,y=mean))+
  #  geom_point(aes(color=`Scoring approach`))+
  geom_line(aes(color=`Scoring approach`))+
  geom_ribbon(aes(ymin=lower,ymax=upper,fill=`Scoring approach`),alpha=0.3)+
  scale_x_log10()+
  labs(y="Sensitivity",x="Maximum rank")+
  theme_bw()+
  ggsave(paste0(save_path_figures,"Sensitivities.pdf"))
