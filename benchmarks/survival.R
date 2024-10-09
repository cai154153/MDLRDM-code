#When using the code in your research work, please cite 
#our paper
#Rappoport, N. and R. Shamir (2018). "Multi-omic and multi-view clustering algorithms: review and cancer benchmark." NUCLEIC ACIDS RESEARCH 46 (20): 10546-10562.

#################################################
get.subtype.survival.path <- function(subtype) {

  datasets.path = get.dataset.dir.path()
  survival.file.path = file.path(datasets.path, subtype, 'survival')
  return(survival.file.path)
}

get.dataset.dir.path <- function() {
  return('Data')
}

########################################################################################################################
#Please set the current file path as working directory, such as setwd('C:/code/benchmarks')
library("survival")
library("ggplot2")
library("ggpubr")
library("survminer")
library(survival)
library(parallel)
library(readxl)
richnum=array();
survivalValue=array();

subtype_clustering <- read_excel("Data/subtypes_aml.xlsx", 
                             col_names = F)
row.names(subtype_clustering)<-subtype_clustering$...1
subtype_clustering<-t(subtype_clustering)
subtype_clustering<-subtype_clustering[-1,]
subtype_clustering<-t(subtype_clustering)
for (i in 1:1)
{
  group=as.integer(subtype_clustering[i,])
  names(group)=colnames(subtype_clustering)
  
}

  groups=group
  subtype="aml"

  survival.file.path = get.subtype.survival.path(subtype)
  uju = read.table(survival.file.path, header = TRUE)
  patient.names = names(groups)
  patient.names.in.file = as.character(uju[, 1])
  patient.names.in.file = toupper(gsub('-', '\\.', substring(patient.names.in.file, 1, 12)))
  
  stopifnot(all(patient.names %in% patient.names.in.file))
  
  indices = match(patient.names, patient.names.in.file)
  mydata = uju[indices,]
  mydata["cluster"] <- group
  mydata$Survival[is.na(mydata$Survival)] = 0
  mydata$Survival=mydata$Survival/365
  mydata$Death[is.na(mydata$Death)] = 0
  fitt<-survfit(Surv(Survival, Death)~cluster, data=mydata);
  #ggsurvplot(fitt,data = mydata,surv.median.line = "hv");
  print(fitt)
  sur_dif<-survdiff(formula=Surv(Survival, Death)~cluster, data=mydata);
  print(sur_dif)
  ggsurvplot(fitt,data=mydata,
             pval = TRUE, 
             conf.int = TRUE,
             legend.title=subtype,#legend.labs=c("Cluster1","Clustehttp://127.0.0.1:24745/graphics/plot_zoom_png?width=951&height=900r2"),
             xlab="Time(Years)",
             risk.table = TRUE, # Add risk table
             risk.table.col = "strata", # Change risk table color by groups
             #surv.median.line = 'hv',
             pval.method = T,
             linetype = "strata", # Change line type by groups
             surv.median.line = "hv", # Specify median survival
             #ggtheme = theme_bw(), # Change ggplot theme
             #palette = c("#E7B800", "#2E9FDF")
  )
  
  # pvalue=get.logrank.pvalue(sur_dif$res)
  survivalValue= -log10(sur_dif$pvalue)
# }
    survivalValue
    
    