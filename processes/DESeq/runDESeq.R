# Run DESeq to find the ramp and peak differential expression for each species


library(DESeq2)
library(stringr)


# set some useful constants:
spcs <- c("BrDi","HoVu","MeNu1","StLa","NaSt")
spcs <- setNames(spcs,spcs)

TtoX <- c(`T-1`=0,T0=1,T1=2,T3=3,T4=4)
TtoF <- c(`T-1`="ramp0",T0="peak0",T1="peak1",T3="ramp1",T4="ramp1")

# read the full expression count tables:

pipelineOutDir <- "/mnt/NOBACKUP/mariansc/share"
tblDir <- file.path(pipelineOutDir,"orthos/exprTbls")

files <- dir(tblDir,pattern=paste0("(",paste(spcs,collapse="|"),")_expected_countTbl.txt"),full.names = T)
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
  
  # contrast tests
  if( "Tframp0"  %in% resultsNames(dds)){
    resRamp <- results(dds, contrast=c("Tf","ramp1","ramp0"))
  } else { # T-1 (ramp0) is missing in StLa... use T0 instead
    resRamp <- results(dds, contrast=c("Tf","ramp1","peak0"))    
  }
  resPeak <- results(dds, contrast=c("Tf","peak1","peak0"))
  
  return(list(vst=m,resRamp=resRamp,resPeak=resPeak))
}) -> DE

save(DE,file="DE.RData")
