library(readr)

# load fasta index
loadFaidx <- function(faidxFile){
  readr::read_tsv(faidxFile, col_names = c("seqID","seqLen","filePos","lineBases","lineBytes"))
}

# convert a group list to a character vector with full sequence names
grpToChar <- function(grp){
  seqsPerSpecies <- sapply(grp,length)
  paste(rep(names(seqsPerSpecies),seqsPerSpecies),unlist(grp),sep="|")
}

# load a set of sequences from fasta file
loadIndexedFasta <- function(seqIDs,fastaFileName,faidx){
  
  fastaFile <- file(fastaFileName,"r")
  
  idx <- match(seqIDs,faidx$seqID)
  
  lapply( setNames(idx,seqIDs), function(i){
    pos <- seek(con = fastaFile,where = faidx$filePos[i])
    cat(pos,"\n")
    paste(readLines(fastaFile,n = ceiling(faidx$seqLen[i]/faidx$lineBases[i])), collapse="")
  }) -> seq
  close(fastaFile)
  
  return(seq)
}


writeFastaFile <- function(outFileName,seq){
  writeLines( paste(paste0(">",names(seq)),seq,sep="\n"), con = outFileName)
}


#### Example:
# 
# faidx <- loadFaidx("/mnt/NOBACKUP/mariansc/share/orthos/orthoMCL/allProteomes.fasta.fai")
# fastaFileName="/mnt/NOBACKUP/mariansc/share/orthos/orthoMCL/allProteomes.fasta"
# 
# # load groups
# source("/mnt/users/lagr/networkSimilarity/R/loadOrthoGroups.R")
# pipelineOutDir <- "/mnt/NOBACKUP/mariansc/share"
# grps <- loadOrthoGrpsArray(orthoGrpFile = file.path(pipelineOutDir,"orthos/orthoMCL/groups.txt"))
# 
# 
# seqIDs <- grpToChar(grps["grp102255",])
# seq <- loadIndexedFasta(seqIDs,fastaFileName,faidx)
# writeFastaFile(outFileName ="~/temp/grp102255.fa", seq)
# 

