rm(list = ls())


load("diet_FDRnotsel.Rdata")



options(max.print=1000000)
diet_data<-diet_data[,-grep("batch36",colnames(diet_data))]


#del<-function(x){sum(is.na(x))/ncol(diet_data)*100} 
#apply(diet_data,1,del)
#diet_data<-diet_data[apply(diet_data,1,del)<80,]

save(diet_data,file="diet_FDRnotsel_without36.Rdata")



load("C:/Users/Cycy/Desktop/diet_R/sample_batch_list_remove36.Rdata")
clinical$diettime=paste0(clinical$diet,clinical$time)
sample_naguide_2=clinical[,c(1,9)]
write.csv(sample_naguide_2,file="sample_naguide_2.csv")


diet_data2<-read.csv("diet_NAguide_data/diet_naguide_result_2.csv",row.names=1)
new_diet_data2 <- diet_data2[order(row.names(diet_data2)), order(colnames(diet_data2))]



diet_data3<-read.csv("diet_NAguide_data/diet_naguide_result_3.csv",row.names=1)
new_diet_data3 <- diet_data3[order(row.names(diet_data3)), order(colnames(diet_data3))]



load("diet_FDRnotsel_remove36_NAguide_NA80_sek.Rdata")
new_diet_data <- diet_data[order(row.names(diet_data)), order(colnames(diet_data))]


all(new_diet_data==new_diet_data3)

library(daff)
d=diff_data(new_diet_data,new_diet_data3)

