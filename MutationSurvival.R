options(stringsAsFactors = F)
library(tidyverse)
library(survival)
MutationSurvial <- function(MutDat, SurDat, gene, Amino_Acid_Change = NULL){
    Sample <- unique(MutDat[[1]])
    SampleCount <- length(Sample)
    if(is.null(Amino_Acid_Change)){
        MutateSample <- MutDat[[1]][MutDat[["gene"]] == gene] 
    }else{
        MutateSample <- MutDat[[1]][MutDat[["gene"]] == gene & 
                                        MutDat[["Amino_Acid_Change"]] == Amino_Acid_Change] 
    }
    MutGroup <- ifelse(Sample %in% MutateSample, 1, 0)
    MutDat1 <- data.frame(sample = Sample, MutGroup = MutGroup)
    JoinDat <- inner_join(MutDat1, SurDat[c(1,2,4)], by = "sample")
    MutCount <- length(MutateSample)
    NonCount <- SampleCount - MutCount
    DeathMuta <- sum(JoinDat$MutGroup == 1 & JoinDat$X_EVENT == 1)
    DeadCount <- sum(JoinDat$X_EVENT == 1)
    DeathMutaPer <- DeathMuta/DeadCount*100
    DeathMutaOutput <- paste0(DeathMuta, "(", round(DeathMutaPer,2), "%)")
    DeathNon <-  sum(JoinDat$MutGroup == 0 & JoinDat$X_EVENT == 1)
    DeathNonPer <- DeathNon/DeadCount*100
    DeathNonOutput <- paste0(DeathNon, "(", round(DeathNonPer,2), "%)")
    if(length(MutateSample) == 0){
        result <- c(MutCount = MutCount, NonCount = NonCount, 
                    DeathMuta = DeathMutaOutput, DeathNon = DeathNonOutput, 
                    HR = NA, CI = NA,
                    p.val = NA)
    }else{
        library(survival)
        my.surv <- Surv(JoinDat$X_TIME_TO_EVENT, JoinDat$X_EVENT)
        ##计算HR以及95%CI
        ##修改分组参照
        data.survdiff <- survdiff(my.surv ~ JoinDat$MutGroup)
        p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
        HR = (data.survdiff$obs[2]/data.survdiff$exp[2])/(data.survdiff$obs[1]/data.survdiff$exp[1])
        up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
        low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
        result <- c(MutCount = MutCount, NonCount = NonCount, 
                    DeathMuta = DeathMutaOutput, DeathNon = DeathNonOutput, 
                    HR = round(HR, 3), CI = paste0(round(low95, 3), "-", round(up95, 3)),
                    p.val = round(p.val, 3))
    }
    return(result)
}