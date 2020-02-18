> 基于UCSC XENA下载的TCGA数据集，批量分析某一个情况的结果

# MutationSurvial

计算某一基因的突变或某一基因的具体突变类型和预后的关系。

## 输入

>MutationSurvial(MutDat, SurDat, gene, Amino_Acid_Change = NULL)

- MutDat: UCSC XENA下载的突变数据

![image-20200122201829098](https://tva1.sinaimg.cn/large/006tNbRwgy1gb5mxopj1pj315d05amzk.jpg)

- SurDat: UCSC XENA下载的生存数据

![image-20200122205931193](https://tva1.sinaimg.cn/large/006tNbRwgy1gb5mzwftctj30m5056gml.jpg)

- gene: 想要检索的基因
- Amino_Acid_Change：基因突变的具体形式。默认是不输入。可以具体制定

## 输出

结果包括7列：

- MutCount： 突变的样本数
- NonCount:   没有突变的样本数
- DeathMuta: 突变当中死亡的样本数
- DeathNon:  没有突变的样本当中死亡的个数
- HR: 使用log-Rank分析预后的HR值(如果突变的没有死亡则为NA)
- CI: 使用log-Rank分析预后的95%CI(如果突变的没有死亡则为NA)
- p.val: 使用log-Rank分析预后的P值(如果突变的没有死亡则为NA)

![image-20200122211855159](https://tva1.sinaimg.cn/large/006tNbRwgy1gb5njtytcoj30o701p0su.jpg)

## 例子

批量计算TCGA当中`p53`的`p.R273H`和预后的关系

```R
# identify the mutation and survival index
mutateIndex <- Sys.glob("~/project/TCGADat/MutateDat/*")
SurvivalIndex <- Sys.glob("~/project/TCGADat/survivalDat/*")
# write a loop to calculate the mutation data
R175H <- c()
for(i in mutateIndex){
    Dat <- read.delim(gzfile(i))
    Res1 <- MutateCount(data = Dat, "TP53", "p.R175H")
    Res1 <- c(CacnerType = str_extract(i, "(?<=-).+?(?=\\.)"),
              Res1)
    R175H <- rbind(R175H, Res1)
}
```

# MutationSurvialPlot

绘制突变和预后相关的生存曲线。生存曲线的颜色采用**默认**的颜色配置。如果需要调整则需要修改代码

## 输入

> MutationSurvialPlot(MutDat, SurDat, gene, 
>                                 Amino_Acid_Change = NULL)

- MutDat(突变数据集): 数据类型**同上**
- SurDat(生存数据集): 数据类型**同上**
- gene: 想要检索的基因
- Amino_Acid_Change: 具体基因的突变

## 输出

![image-20200217193808304](https://tva1.sinaimg.cn/large/0082zybply1gbzmsg194rj30nd0n2jsm.jpg)

## 例子：

绘制👆有突变的肿瘤的生存曲线

```R
## based on the above result, we plot the survivial curve wihch cancer has the mutation sample
R273HPlot <- list()
for(i in rownames(R273H1)[!is.na(R273H1[,5])]){
    Mutindex <- str_subset(mutateIndex, i)
    Surindex <- str_subset(SurvivalIndex, i)
    MutDat <- read.delim(gzfile(Mutindex))
    SurDat <- read.delim(gzfile(Surindex))
    R273HPlot[[i]] <- MutationSurvialPlot(MutDat = MutDat, SurDat = SurDat, gene = "TP53", 
                              Amino_Acid_Change = "p.R273H")
    ggsave(paste0(i, "_R273H.pdf"), width = 4, height = 4)
}
```



# ExprExtract

提取某一些基因的表达量，同时基于样本好区分癌和正常。(这里把所有以0..结尾的都当作了癌症)

## 输入

>ExprExtract(ENSID, cancerid, index)

- ENSID: 这个是基于genecode22的基因ID号，可以是**多个编号**。例如: ENSG00000096088.15

- cancerid: pan-cancer当中目标肿瘤的简写。

- index: 制定数据存放的位置

  **需要注意的是：**

- tpmindex和canerid之间存在一定的包含的关系。例如tpmindex是：/Users/lihao/project/TCGADat/tpmDat/ACC_tpm.csv.gz。而cancerid则是: ACC。我们是通过cancerid查找唯一的TPMindex。

- tpm的表达数据是csv的压缩文件。

## 输出

结果包括至少三列：

- 目标基因的表达量

- 肿瘤类型

- 肿瘤的分组: cancer OR normal 

  ![image-20200216231857568](https://tva1.sinaimg.cn/large/0082zybpgy1gbynih2eouj30jn055myr.jpg)

## 例子

提取基因`ENSG00000096088.15`的表达量:

```R
tpmindex <- Sys.glob("~/project/TCGADat/tpmDat/*")
cancerid <- str_extract(tpmindex, "(?<=tpmDat\\/).+?(?=_)")
PGCExpr <- GeneExprPan("ENSG00000096088.15", cancerid)
```

# PanExprSummary

summary the expression data extracted from `ExprExtract` function. this is a script. not a function.

the summary contains a table including the `mean, mdian and sd` in every cancer and a box plot group by sampe type.