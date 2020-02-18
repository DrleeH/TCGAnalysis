### 需要提前制定三个变量
### 1. 需要提取的基因的ENSID(genecode22):ENSG00000096088.15
### 2. 所有pancancer表达数据的存储位置
### 3. 所有cancer的简写(和存储文件命有一定的惯性)
options(stringsAsFactors = F)
library(tidyverse)
ExprExtract <- function(ENSID, cancerid, index){
    dat <- c()
    for(i in cancerid){
        tpmindex1 <- str_subset(index, i)
        tpmdat <- vroom::vroom(gzfile(index))
        tpmdat <- tpmdat %>% as.data.frame() %>% column_to_rownames(var = "...1")
        TargetDat <- t(tpmdat[ENSID, ]) %>% as.data.frame()
        TargetDat$CacnerType <- i
        TargetDat$Group <- ifelse(str_detect(rownames(TargetDat), "0..$"), "Cancer", "Normal")
        dat <- rbind(dat,TargetDat)
        
    }
    return(dat)
}