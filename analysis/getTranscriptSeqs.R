####
#
# get the sequence of the "conserved" genes
#


#
# get grpIDs for the genes that are significantly DE in all species
#

source("~/GRewd/pipeline/R/orthoGrpTools.R")
DEmat <- readRDS("~/GRewd/pipeline/data/DEmat.RDS")  # from: superGeneModel.Rmd
dim(DEmat$peak$pVal)

sigMatUp <- lapply(DEmat, with, { ifelse(is.na(pAdj), F, pAdj < 0.05 & lfc > 1) })
sigMatDown <- lapply(DEmat, with, { ifelse(is.na(pAdj), F, pAdj < 0.05 & lfc < -1) })

# get genes with significant peak/ramp up/down in all species
grpIDs <- unique(c(names(which(apply(sigMatUp$ramp,1,all))),
                   names(which(apply(sigMatUp$peak,1,all))),
                   names(which(apply(sigMatDown$ramp,1,all))),
                   names(which(apply(sigMatDown$peak,1,all)))))


#
# get sequences
#


####
#
# extractFromFasta
#
library("RLinuxModules")
moduleInit()
module("load samtools")

extractFromFasta <- function(fastaFile, seqIDs, outFile=NULL){
  cmd <- paste("samtools faidx", fastaFile, paste(shQuote(seqIDs),collapse=" "))
  if(is.null(outFile)){
    return(system(cmd,intern = T))
  } else {
    system(paste(cmd,">",outFile))
  }
}

####

# get seqIDs from orthoGrps
orthoPath <- "/mnt/NOBACKUP/mariansc/share/orthos"
grps <- loadOrthoGrpsArray(file.path(orthoPath,"splitGroups/goodGroups.txt"))



spcs <- c("BrDi","HoVu","MeNu1","NaSt","StLa")
spcs <- setNames(spcs,spcs)

# look up the table to get the original transcriptIDs
seqID2transcriptID <- function(seqIDs, spc){
  filename <- file.path(orthoPath,"longestORFs",paste0(spc,".tbl"))
  seqIDtbl <- sapply(seqIDs,function(seqID){
    system(paste("grep",seqID,filename),intern = T)
  })
  sapply(strsplit(seqIDtbl,split = "\t"), function(x){
    sub("cds\\.(TR[0-9]+\\|c[0-9]+_g[0-9]+_i[0-9]+)\\|m\\.[0-9]+","\\1",x[2])
  })
}


#
# Store sequence
#
outDir = "seqs"
dir.create(outDir)

# get alignment files:
alnFiles <- file.path(orthoPath,"grpAligned",paste0(sub("\\.[0-9]+","",grpIDs),".aln"))
cdsalnFiles <- file.path(orthoPath,"pal2nal",paste0(sub("\\.[0-9]+","",grpIDs),".cds.aln"))

# copy them to the outDir
file.copy(alnFiles,outDir)
file.copy(cdsalnFiles,outDir)


# for each grp
lapply(setNames(grpIDs,grpIDs),function(grpID){
  # create path for grp
  grpDir <- file.path(outDir,grpID)
  dir.create(grpDir)
  
  #for each spc in grp
  lapply(spcs,function(spc){
    seqIDs <- grps[[grpID,spc]]
    transcriptIDs <- seqID2transcriptID(seqIDs,spc)
    fastaFile <- file.path("/mnt/NOBACKUP/mariansc/share/trinity",spc,paste0(spc,".fasta"))
    
    # for each paralogous seq
    for(i in seq_along(transcriptIDs)){
      outFile <- file.path(grpDir,paste0(spc,"_",names(transcriptIDs)[i],".fa"))
      extractFromFasta(fastaFile = fastaFile, seqIDs = transcriptIDs[i],
                       outFile = outFile)
    }
  })
})

