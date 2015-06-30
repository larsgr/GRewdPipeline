bresFileName="/mnt/NOBACKUP/mariansc/share/orthos/orthoMCL/all_vs_all.out"
# outFile="all_vs_all.out.index.RData"

indexBlastResults <- function(bresFileName,outFile){
  fileSize <- file.info(bresFileName)$size

  bresFile <- file(bresFileName,"r")
  qseq.prev <- "FOOBAR"
  filepos <- 0
  n=100000 # number of lines to read per iteration
  idx.qseq <- character(0)
  idx.pos <- numeric(0)
  repeat{
    x <- readLines(con = bresFile,n = n)
    if(length(x)==0) break
    
    qseq <- str_extract(x,"^[^\t]+")
    seqPos <- cumsum(c(filepos,nchar(x)+1))
    newSeqsIdx <- which(qseq != c(qseq.prev, qseq[1:(length(x)-1)]))
    
    idx.qseq <- c(idx.qseq, qseq[newSeqsIdx])
    idx.pos <- c(idx.pos, seqPos[newSeqsIdx])
    
    qseq.prev <- qseq[length(x)]
    filepos <- seek(con = bresFile)
    cat(floor(100*filepos/fileSize),"%","\r",sep = "")
  }
  
  close(bresFile)
  
  cat("\nWriting index to file:",outFile)
  bresIdx <- data.frame(qseq=idx.qseq,pos=idx.pos,stringsAsFactors=F)
  save(bresIdx,file = outFile)
}



loadIndexedBres <- function(seqId,bresFileName,bresIdx){
  
  bresFile <- file(bresFileName,"r")
  
  idx <- match(seqId,bresIdx$qseq)
  
  lapply(setNames(idx,seqId), function(i){
    seek(con = bresFile,where = bresIdx$pos[i])
    len <- bresIdx$pos[i+1] - bresIdx$pos[i]
    # TODO: handle i+1 > length
    x <- readChar(con = bresFile, nchars = len)
    
    read.table(text = x, sep="\t", stringsAsFactors=F,
               col.names = c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", 
                             "qstart", "qend", "sstart", "send", "evalue", "bitscore"))
  }) -> bres
  close(bresFile)
  
  return(bres)
}


# seqId = sample(bresIdx$qseq,100)
# 
# system.time( bres <-loadIndexedBres(seqId,bresFileName,bresIdx) )


