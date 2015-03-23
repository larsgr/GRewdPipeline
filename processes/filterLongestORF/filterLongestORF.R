library(seqinr)
library(data.table)


filterFastaLongestORF <- function(inFasta, outFasta, outTable){

  ORF <- read.fasta(inFasta)
  
  # create table with sequence length, name and name with isoform removed
  dt <- data.table(len = getLength(ORF), name = names(ORF),
                   nameNoIso = unlist(lapply(strsplit(names(ORF),"_i"),"[",1)))

  # For each "gene", get rank order of the length of its isoforms
  dt[,isoLengthRank := rank(-len,ties.method = "first"),by=nameNoIso]


  # Get names of the longest isoforms
  longestIso <- dt$name[ dt$isoLengthRank==1 ]
  
  # change names e.g "cds.TR100006|c0_g1_i1|m.887688" -> TR100006_c0_g1
  longestIsoNewName <- gsub("\\|","_",substring(dt$nameNoIso[ dt$isoLengthRank==1 ],5))
  
  write.table(cbind(longestIsoNewName,longestIso), file = outTable,
              quote = F,row.names = F, col.names = F, sep="\t")
  
  # Write fasta sequence of the longest isoforms
  write.fasta(sequences=ORF[longestIso], names=longestIsoNewName, file.out=outFasta)
}

args <- commandArgs(trailingOnly = T)
if( length(args) != 3){
  cat("Error! Invalid number of parameters!\n\n")
  cat("Usage: Rscript filterLongestORF.R <inFasta> <outFasta> <outTable>\n")
  cat("
The <inFasta> file must have sequence names in the format created by running transdecoder
on a trinity assembly e.g.: cds.TR100006|c0_g1_i1|m.887688
The <outFasta> file will only contain the longest isoform of each gene and will be
given a new simplified name like: TR100006_c0_g1
The <outTable> file contains the mapping between the new names and original names.
")
} else {
  filterFastaLongestORF( args[1], args[2], args[3] )
}
