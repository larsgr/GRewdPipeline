# source("http://bioconductor.org/biocLite.R")
# biocLite("DESeq2")

library(DESeq2)


pipelineOutDir <- "/mnt/NOBACKUP/mariansc/share"
grps <- loadOrthoGrpsArray(orthoGrpFile = file.path(pipelineOutDir,"orthos/orthoMCL/groups.txt"))

tblDir <- file.path(pipelineOutDir,"orthos/exprTbls")

grp11Cnt <- read.table(file = file.path(tblDir,"grpSingletonCountTbl.txt"), 
                       stringsAsFactors = F)
grpSumCnt <- read.table(file = file.path(tblDir,"grpSumCountTbl.txt"), 
                       stringsAsFactors = F)


counts <- round(grp11Cnt[,!grepl("^wc_",colnames(grp11Cnt))])

# define metaData
pDat <- as.data.frame(str_split_fixed(colnames(counts),"\\.",4)[ ,1:3])
names(pDat) <- c("spcPop","timePoint","mix")
pDat$spcAsm <- as.factor(str_extract(colnames(counts),paste(spcs,collapse="|")))


# get data only for one species

HoVuCnt <- read.table(file = file.path(tblDir,"HoVu_expected_countTbl.txt"), 
                      stringsAsFactors = F)

countData <- round(HoVuCnt[,!grepl("^wc_",colnames(HoVuCnt))])

colData <- as.data.frame(str_split_fixed(colnames(countData),"\\.",4)[ ,1:3])
#colData <- as.data.frame(str_split_fixed(sample(colnames(countData)),"\\.",4)[ ,1:3])
names(colData) <- c("spcPop","timePoint","mix")
colData$time <- TtoX[ colData$timePoint] - 1
colData$peak <- colData$time==1

dds <- DESeqDataSetFromMatrix(countData = countData,
                              colData = colData,
                              design = formula(~ peak + time))

# estimate size factors
dds <- estimateSizeFactors(dds)

cor.test(sizeFactors(dds),colData$time) # corralates with time?

# estimate dispersion
dds <- estimateDispersions(dds)
plotDispEsts(dds)


ddsPeak <- nbinomLRT(dds, reduced = ~ time)

ddsTime <- nbinomLRT(dds, reduced = ~ peak)

idx <- match(gsub("_i.+$","",gsub("\\|","_",Coldbres)),rownames(countData))
len <- nrow(countData)

plot(rank(results(ddsPeak)$stat)[idx]/len,rank(results(ddsTime)$stat)[idx]/len, pch=20,
     xlim=c(0,1),ylim=c(0,1), 
     col=ifelse(pmin(results(ddsPeak)$padj[idx],results(ddsTime)$padj[idx])<0.05,"red","black"))

plot(rank(results(ddsPeak)$stat),rank(results(ddsTime)$stat), pch=20,
     col=ifelse(rownames(countData) %in% Cold.grp,"red","black"))

hist(log10(results(ddsTime)$stat))

length(results(ddsTime)$pvalue)

sum(results(ddsPeak)$pvalue<0.00001,na.rm = T)/sum(!is.na(results(ddsPeak)$pvalue))
sum(results(ddsPeak)$pvalue>0.001,na.rm = T)

plot(results(ddsTime)$pvalue,results(ddsPeak)$pvalue,log="xy",xlim=c(1e-10,1),ylim=c(1e-10,1))

res[order(res$padj), ]
res <- results(ddsTime)
i=1
i=i+5
plotTExp(log2(1+countData[rownames(res[order(res$padj), ])[i:(i+1)],]))

dds <- DESeq(dds,reduced = ~ time)
dds <- DESeq(dds,reduced = ~ peak)
resPeak <- results(DESeq(dds,reduced = ~ peak))
resTime <- results(DESeq(dds,reduced = ~ time))

identical(resPeak,resTime)
resPeak[order(resPeak$padj), ]
resTime[order(resTime$padj), ]

plotTExp(log2(1+countData[rownames(resPeak[order(resPeak$padj), ])[1:5],]))
plotTExp(log2(1+countData[rownames(resTime[order(resTime$padj), ])[1:5],]))

plotMA(res)
?nbinomLRT
