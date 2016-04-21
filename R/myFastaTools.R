# Read fasta file as named character vector
readFasta <- function(inFastaFile){
  # read fasta file
  txt.fasta <- readLines(inFastaFile)
  # convert to into named character vector
  headerLines <- grepl("^>",txt.fasta)
  seqIDs <- sub("^>","",txt.fasta[headerLines]) # OBS assumes only seqID on header
  seqs <- tapply(txt.fasta[!headerLines],
                 seqIDs[cumsum(headerLines)[!headerLines]],
                 paste,collapse="")
  return(seqs)
}

# Write named character vector as fasta file
writeFasta <- function(seqs,outFastaFile){
  write(paste0(">",names(seqs),"\n",seqs,collapse = "\n"), 
        file = outFastaFile )
}

# # remove columns with gaps in all sequences
# removeGapInAll <- function(seqs){
#   seqArray <- strsplit(seqs, split="")
#   gapMat <- sapply(seqArray, "==","-")
#   gapColumns <- apply(gapMat,1,all)
#   lapply(seqArray, function(seq){paste(seq[!gapColumns],collapse="")})
# }