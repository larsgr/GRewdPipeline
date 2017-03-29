# Run DESeq to generate variance stabilized expression for all


library(DESeq2)
library(stringr)


# set some usefull constants:
spcs <- c("BrDi","HoVu","MeNu1","MeNu2","StLa","NaSt")
spcs <- setNames(spcs,spcs)

TtoX <- c(`T-1`=0,T0=1,T1=2,T3=3,T4=4)
TtoF <- c(`T-1`="ramp0",T0="peak0",T1="peak1",T3="ramp1",T4="ramp1")

# read the full expression count tables:

pipelineOutDir <- "/mnt/NOBACKUP/mariansc/share"
tblDir <- file.path(pipelineOutDir,"orthos/exprTbls")

files <- dir(tblDir,pattern="full_expected_countTbl.txt",full.names = T)
exprCnt <- lapply(setNames(files,str_extract(basename(files),"^[^_]+")),
                  read.table, stringsAsFactors = F)





# for each species:
lapply(exprCnt,function(countData){

  # fix the T-1 name
  names(countData) <- sub("T\\.1","T-1",names(countData))
  
  colData <- as.data.frame(str_split_fixed(names(countData),"\\.",4)[ ,1:3])
  names(colData) <- c("spcPop","timePoint","mix")
  colData$Tf <- as.factor(TtoF[as.character(colData$timePoint)])
  
  
  dds <- DESeqDataSetFromMatrix(countData = round(countData),
                                colData = colData,
                                design = formula(~ Tf ))
  
  dds <- DESeq(dds)
  
  vsd <- varianceStabilizingTransformation(dds)
  m <- assay(vsd)
  colnames(m) <- names(countData)
  
  return(m)
}) -> VST

save(VST,file="VST.RData")
