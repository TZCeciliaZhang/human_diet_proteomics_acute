rm(list = ls())

library("readxl")
library("tidyverse")
library("stringr")
library("vioplot")
library("openxlsx") 
library("dplyr")



load("batch_list.Rdata")

result=list()
for (i in 1:length(list)) {
  list[[i]]%>%
  
    select(contains("Accession")|contains("Protein FDR Confidence: Combined")|contains("Abundance Ratio: "))%>%
    as.data.frame()->batch
  rownames(batch)=batch[,1]
  batch<-batch[,-1]

  FDRsel<-which(batch$`Protein FDR Confidence: Combined`!="Low")
  batch<-batch[FDRsel,]%>%
    select(contains("Accession")|contains("Abundance Ratio: "))->batch
 
    type_i<-str_replace(list.files("diet_protein_alldata/",pattern = ".xlsx")[i],".xlsx","")

    colnames(batch)<-c(paste0(type_i,str_replace(colnames(batch),"Abundance Ratio: ","AR"),sep=""))
 
    Accession=rownames(batch)
    batch=cbind(Accession,batch)
 
    result_name=paste0(type_i,"FDRsel.xlsx",sep="")
    result[[i]]<-write.xlsx(batch,result_name,rowNames=F,colNames=T) }


list.files("diet_protein_alldata_FDRsel/",pattern = ".xlsx",full.names = TRUE)%>%
  lapply(readxl::read_excel)%>%
  reduce(full_join,by="Accession") -> diet_data

options(max.print=1000000)
diet_data<-as.data.frame(diet_data)
rownames(diet_data)=diet_data[,1]
diet_data<-diet_data[,-1]

save(diet_data,file="diet_data.Rdata")
write.xlsx(diet_data,"diet_data.xlsx",rowNames=T,colNames=T)



load("batch_list.Rdata")

result=list()
for (i in 1:length(list)) {
  list[[i]]%>%

    select(contains("Accession")|contains("Abundance Ratio: "))%>%
    as.data.frame()->batch
  rownames(batch)=batch[,1]
  batch<-batch[,-1]

  type_i<-str_replace(list.files("diet_protein_alldata/",pattern = ".xlsx")[i],".xlsx","")
 
  colnames(batch)<-c(paste0(type_i,str_replace(colnames(batch),"Abundance Ratio: ","AR"),sep=""))

  Accession=rownames(batch)
  batch=cbind(Accession,batch)

  result_name=paste0(type_i,"FDRnotsel.xlsx",sep="")
  result[[i]]<-write.xlsx(batch,result_name,rowNames=F,colNames=T) }


list.files("diet_protein_alldata_FDRnotsel/",pattern = ".xlsx",full.names = TRUE)%>%
  lapply(readxl::read_excel)%>%
  reduce(full_join,by="Accession") -> diet_data

options(max.print=1000000)
diet_data<-as.data.frame(diet_data)
rownames(diet_data)=diet_data[,1]
diet_data<-diet_data[,-1]


colnames(diet_data)<-gsub('[(]', '.', colnames(diet_data))
colnames(diet_data)<-gsub('[)/]', '', colnames(diet_data))
colnames(diet_data)<-gsub(' ', '', colnames(diet_data))


options(max.print=1000000)
diet_data<-diet_data[-grep("Q5T0Z8",rownames(diet_data)),]


save(diet_data,file="diet_FDRnotsel.Rdata")






  





