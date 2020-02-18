> Based on TCGA datasets download  from UCSC XENA. Some useful functions were writted based on R to do the Cancer related analysis.

# MutationSurvial

calculate the correaltion between one gene mutation or specific mutation type of a gene and survival

## import

> MutationSurvial(MutDat, SurDat, gene, Amino_Acid_Change = NULL)

- MutDat: mutation data download from UCSC XENA 

![image-20200122201829098](https://tva1.sinaimg.cn/large/0082zybply1gbzm5joeppj315d05ata8.jpg)

- SurDat: Survival information download from UCSC XENA

![image-20200122205931193](https://tva1.sinaimg.cn/large/0082zybply1gbzm689gugj30m5056mxr.jpg)

- gene: gene which we wanted searched 
- Amino_Acid_Change: the specific mutation. NULL is default.

## output

the result contains 7 column

- MutCount: mutation sample number
- NonCount: non-mutation sample number
- DeatMuta: Death number in mutation group
- DeathNon: Death number in non-mutaion group
- HR: HR value in the survival analysis
- CI: 95% CI value in the survival anlysis
- p.val: p value in the survial analysis

**PS:** if the the sample has no mutaion the result is NA.

![image-20200122211855159](https://tva1.sinaimg.cn/large/0082zybply1gbzmd1n20gj30o701pmx4.jpg)

## example

Survival analysis the mutation type  ``p.R273H`` in  `TP53` gene 

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

plot the survival curve for one gene mutaion or specific mutaion type of a gene. the plot function use the surviminer::ggsurvivalplot

## Import

> MutationSurvialPlot(MutDat, SurDat, gene, 
>                                 Amino_Acid_Change = NULL)

- MutDat(mutation datasets):  the data format is the same above
- SurDat(survival datasets):  the data format is the same above
- gene: one gene we wanted searched
- Amino_Acid_Change: one mutation type of the gene

## Export

![image-20200217193808304](https://tva1.sinaimg.cn/large/0082zybply1gbzmqz6rj7j30nd0n2myg.jpg)

## Example:

plot the survival curve for the result above

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

extract one gene tpm/fpkm/count expression data, making out the sample type(cancer or normal) at the same time.

PS: base on the TCGA sample id, we regared the `0..$` sample as cancer sample

## Import

>ExprExtract(ENSID, cancerid, index)

- ENSID: gene ensemble ID based on genecode 22. we can import **multiple** id numbers. 例如: ENSG00000096088.15

- cancerid:  cancer type abbrevition

- index: expresssion data stored location

  ![image-20200217195826260](https://tva1.sinaimg.cn/large/0082zybply1gbznc3jk2ij30h4012wel.jpg)

  **PS**

- expression data should be the gz compressed file.

## Export

Result contains as least three columns：

- Expression of candidate genes

- Cancer type

- sample group( cancer OR normal )

  ![image-20200216231857568](https://tva1.sinaimg.cn/large/0082zybply1gbzn25belqj30jn055jry.jpg)

## Example

提取基因`ENSG00000096088.15`的表达量:

```R
tpmindex <- Sys.glob("~/project/TCGADat/tpmDat/*")
cancerid <- str_extract(tpmindex, "(?<=tpmDat\\/).+?(?=_)")
PGCExpr <- GeneExprPan("ENSG00000096088.15", cancerid, tpmindex)
```

# PanExprSummary

summary the expression data extracted from `ExprExtract` function. this is a script. not a function.

the summary contains a table including the `mean, mdian and sd` in every cancer and a box plot group by sampe type.