setwd('D:/DrugMatrix/practice/')
library(dplyr)
library(sqldf)
library(data.table)
library("annotate")
library("rat2302.db")

compoundAnnotation <- read.csv.sql('./COMPOUND_ANNOTATIONS.txt',sql = "select CAS_NUM, NAME from file", sep="\t")
affyLiver <- read.table('./AFFY_LIVER.txt',header=TRUE,check.names = FALSE)
affyAnnotation <- read.csv.sql('./Affymetrix_annotation.txt',sql="select EXPERIMENT, EXPERIMENT_NAME, TISSUE_NAME, COMPOUND_NAME, DOSE, DOSE_UNIT, ADMINISTRATION_ROUTE,ADMINISTRATION_FREQUENCY, TIME, TIME_UNIT from file",sep="\t")
#affyLiverAnnotation<-filter(affyAnnotation,TISSUE_NAME=="LIVER")
affyAnnotation<-affyAnnotation %>% select(c(EXPERIMENT,EXPERIMENT_NAME,COMPOUND_NAME))

affyGenes <- read.table('./Affymetrix_genes.txt', sep="\t", header=TRUE,check.names=FALSE)
histoPathology<-read.table('./HISTOPATHOLOGY_EXP_REPORT.txt',sep="\t",header=TRUE)
LiverHisto<-filter(histoPathology,TISSUE_NAME=="LIVER")
LiverHisto<-LiverHisto %>% select(c(EXPERIMENT,AVERAGE_SEVERITY_TOTAL,AVERAGE_SEVERITY_AFFECTED,PVALUE,PERCENT_INCIDENCE))
affyLiverHisto<-merge(LiverHisto,affyAnnotation)
affyLiverHisto<- affyLiverHisto %>% group_by(COMPOUND_NAME)%>%summarize(AVERAGE_SEVERITY_TOTAL=sum(AVERAGE_SEVERITY_TOTAL),AVERAGE_SEVERITY_AFFECTED=sum(AVERAGE_SEVERITY_AFFECTED),PVALUE=min(PVALUE),PERCENT_INCIDENCE=sum(PERCENT_INCIDENCE))

#files<-list.files(path="./test_exp/",recursive=TRUE,full.names=TRUE)
#pul<-read.csv.sql(files[1],sql="select GENE_NAME, ADJUSTED_LOG_RATIO,STDEV_OF_LOG_RATIO,SCORE,INTENSITY,ORGANISM_DESCRIPTION,COMPOUND_NAME from file", sep="\t")
#write.csv(file = 'Affy_probes.txt',affyLiver[,1],row.names = FALSE,quote=FALSE)

gexData <- merge(affyGenes,affyLiver)
#names(gexData)[match(affyAnnotation[,1],names(gexData))]<-affyAnnotation[,2]
drops<-c("PROBE","Species","Gene_Name")
gexData<-gexData[,!(names(gexData) %in% drops)]

cggex<-setNames(data.frame(t(gexData[,-1])),gexData[,1])
cggex[,"EXPERIMENT"]<-c(as.numeric(rownames(cggex)))
cggex<-cggex[,c(20759,1:20758)]
cggex[,"COMPOUND_NAME"]<-affyAnnotation[match(cggex$EXPERIMENT,affyAnnotation$EXPERIMENT),3]
cggex<-cggex[,c(20760,1:20759)]
set(cggex,j="EXPERIMENT",value=NULL)
cggex<-cggex%>%group_by(COMPOUND_NAME)%>%summarize_each(funs(mean))
cggex<-merge(cggex,affyLiverHisto)

write.table(file = './Affy_Liver.tsv',cggex,row.names = TRUE,col.names=NA, quote=FALSE,sep="\t")

#experiments<-names(gexDat)[!(names(gexDat) %in% "UniProt")]
#experiments<-data.frame(EXPERIMENT=c(experiments))
#expAnnotation<-merge(affyAnnotation,experiments,by.x='EXPERIMENT',by.y='Experiment')

write.table(file='Affy_Liver.tsv',gexData,row.names=FALSE,quote=FALSE,sep="\t")