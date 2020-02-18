MutationSurvialPlot <- function(MutDat, SurDat, gene, 
                                Amino_Acid_Change = NULL){
    library(survival)
    library(tidyverse)
    library(survminer)
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
    if(is.null(Amino_Acid_Change)){
        LegendLab <- c(paste0("MutCount(", MutCount, ")"),
                       paste0("NonCount(", NonCount, ")"))
    }else{
        LegendLab <- c(paste0(Amino_Acid_Change,"(", MutCount, ")"),
                       paste0(Amino_Acid_Change,"(", NonCount, ")"))
    }
    my.surv <- Surv(JoinDat$X_TIME_TO_EVENT, JoinDat$X_EVENT)
    fit <- survfit(my.surv ~ JoinDat[[2]])
    ##计算HR以及95%CI
    ##修改分组参照
    data.survdiff <- survdiff(my.surv ~ JoinDat$MutGroup)
    p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
    ggsurvplot(fit, data = JoinDat ,
               #ggtheme = theme_bw(), #想要网格就运行这行
               conf.int = F, #不画置信区间，想画置信区间就把F改成T
               #conf.int.style = "step",#置信区间的类型，还可改为ribbon
               censor = T, #不显示观察值所在的位置
               #palette = c("#D95F02","#1B9E77"), #线的颜色对应高、低
               
               legend.title = gene,#基因名写在图例题目的位置
               font.legend = 11,#图例的字体大小
               #font.title = 12,font.x = 10,font.y = 10,#设置其他字体大小
               
               #在图例上标出高低分界点的表达量，和组内sample数量
               legend.labs=LegendLab,
               
               #在左下角标出pvalue、HR、95% CI
               #太小的p value标为p < 0.001
               pval = round(p.val, 3))
}