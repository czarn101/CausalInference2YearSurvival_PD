library(data.table)
library(nseval)
library(plyr)
library(dplyr)
library(nnet)
library("anytime")
library(MASS)
library(parallel)
library(doParallel)
library(pacman)

lapply2=function(...,core_num =30,p=NULL,e=NULL,isPar2=T,ee=caller(environment())){
  if(!isPar2){
    return(lapply(...))
  }else if(Sys.info()["sysname"] %in% "Windows"){
    cl <- makeCluster(core_num)#getOption("cl.cores",core_num))
    clusterExport(cl,c(e,"p_load"),envir=ee)
    if(!is.null(p)){
      clusterCall(cl, function() do.call(p_load,as.list(p)))
    }
    r = parLapply(cl,...)
    stopCluster(cl)
  }else{
    r=mcapply(...,mc.cores = core_num)
  }
  return(r)
}

L2 = lapply2(1:1000,function(x) x^2 )

filePath = "C:\\Users\\Albert\\Dropbox\\AlbertPierce\\ParkinsonData\\GoodDataset\\ParkinsonSurvivalDataset_2Treatments_2Year.csv"
dre1 <- data.frame(matrix(ncol = 4, nrow = 0))

df = fread(filePath, header = TRUE)
samples = 10000
causal_survival = function(mydf,i){
  df$V1 = 1:nrow(df)
  rs = as.data.frame(sample(1:nrow(df), nrow(df), replace=TRUE))
  colnames(rs)<-"V1"
  df <- merge(x = rs,y = df, by=c("V1"), all.x=TRUE)
  
  
  df = df[,-c(1)]
  df$personid = 1:nrow(df)
  
  
  df = df[df$Treatment == "No Treatment" | df$Treatment == "carbidopa-levodopa",]
  df$i = as.numeric(ifelse(df$Treatment == "No Treatment", 0, 1))
  
  df = dplyr::select(df, c("personid","survived","Gender","MaritalStatus" ,           
                           "Race","AgeAtDiagnosis",              
                           "Diabetes","OverWeightAndObesity","CerebralInfarction" ,      
                           "AcuteKidneyFailure","HypertensiveDiseases","HeartDisease","i"))
  
  base.model = glm(as.factor(i)~., family=binomial, data = df[,-c(1,2)])
  #VGAM PACKAGE, N NET(multinumial logisitc)
  predictdf <- predict(base.model, newdata = df[,-c(1,2)], type = "response")
  
  predictdf <- as.data.frame(predictdf)
  predictdf <- cbind(predictdf, df[,c(1,2,13)])
  #predictdf[is.na(predictdf)] = 0
  
  predictdf$ipw = ifelse(predictdf$i == 0, 1/(1-predictdf$predictdf), 1/(predictdf$predictdf))
  
  predictdf$ey = predictdf$survived*predictdf$ipw
  my_vector = c()
  
  #Horwitz-Thompson IP weighted estimators
  for(ii in 1:2) { 
    
    my_vector[ii] <- sum(as.numeric(predictdf[predictdf$i == ii-1, ]$ey)/nrow(predictdf))
    
  }
  my_vector
  
  df$i <- as.factor(df$i)
  base.model2 = glm(survived ~.,family = binomial, data = df[,-1])
  df1=df
  
  
  df$i <- "0"
  df1$i <- "1"
  
  
  dftempMerged <- do.call("rbind", list(df, df1))
  dftempMerged$p1 <- predict(base.model2,newdata = dftempMerged[,-1],type = "response")
  dftempMerged1 = dplyr::select(dftempMerged, c("personid","i","p1"))
  TdftempMerged1 <- dcast(setDT(dftempMerged1), personid ~i,fun.aggregate =sum, value.var = "p1")
  colnames(TdftempMerged1) = c("personid","T0","T1")
  
  final <- merge(x = predictdf,y = TdftempMerged1, by=c("personid"), all.x=TRUE)
  final$ft = ifelse(final$i == 0, as.numeric(final$ipw*final$T0), NA)
  final$ft = ifelse(final$i == 1, as.numeric(final$ipw*final$T1), final$ft)
  
  final1 = dplyr::select(final, c("personid","i","ft"))
  
  final2 <- dcast(setDT(final1), personid ~i, fun.aggregate =sum,value.var = "ft")
  final2[is.na(final2)] = 0
  
  dre1 = rbind(dre1,c(i,"IP weighted",my_vector[1],my_vector[2]))
  dre1 = rbind(dre1,c(i,"Y for patients",mean(final2$`0`),mean(final2$`1`)))
  dre1 = rbind(dre1,c(i,"standardized estimators",mean(final$T0),mean(final$T1)))
  dre1 = rbind(dre1,c(i,"doubly robust estimators",
                      (my_vector[1]+mean(final$T0)-mean(final2$`0`)),
                      (my_vector[2]+mean(final$T1)-mean(final2$`1`))))
  
  dre1 = rbind(c(i,"IP weighted",my_vector[1],my_vector[2]),c(i,"Y for patients",mean(final2$`0`),mean(final2$`1`)),c(i,"standardized estimators",mean(final$T0),mean(final$T1)),c(i,"doubly robust estimators",(my_vector[1]+mean(final$T0)-mean(final2$`0`)),(my_vector[2]+mean(final$T1)-mean(final2$`1`))))
  return(dre1)
}

L=lapply2(1:10000,function(x) causal_survival(df,x),p=c("data.table","plyr","dplyr","nnet","anytime","MASS"),e=c("df","causal_survival"))
saveRDS(L,"C:\\Users\\Albert\\Dropbox\\AlbertPierce\\ParkinsonData\\L.RDS") 