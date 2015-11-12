library(seqinr)
library(data.table)
library(stringr)


filterFastaLongestORF <- function(inFasta, outFasta, outTable){

  ORF <- read.fasta(inFasta)
  
  # Convert sequence names
  # remove isoform part and pipes
  if( grepl("_i",names(ORF)[1] )){
    # These are our transcriptomes
    # change names e.g "cds.TR100006|c0_g1_i1|m.887688" -> "TR100006_c0_g1"
    nameNoIso <- sub("cds.(TR[0-9]+)\\|(c[0-9]+_g[0-9]+)_i.*","\\1_\\2",names(ORF),perl="T")
  } else if( grepl("^cds\\.ENA",names(ORF)[1])) {
    # this is the LoPe transcriptome
    # format e.g.: cds.ENA|GAYX01000005|GAYX01000005.1|m.9
    nameNoIso <- sub("cds.ENA\\|(GAYX[0-9]+)\\|.*","\\1",names(ORF),perl="T")
  } else {
    # this is the lolium falster transcriptome
    # change names e.g "cds.f19713_c0s2|m.2" -> "f19713_c0s2"
    nameNoIso <- sub("cds.(f[0-9]+_c[0-9]+s[0-9]+)\\|.*","\\1",names(ORF),perl="T")
  }
  
  # create table with sequence length, name and name with isoform removed
  dt <- data.table(len = getLength(ORF), name = names(ORF),
                   nameNoIso = nameNoIso)

  # For each "gene", get rank order of the length of its isoforms
  dt[,isoLengthRank := rank(-len,ties.method = "first"),by=nameNoIso]


  # Get names of the longest isoforms
  longestIso <- dt$name[ dt$isoLengthRank==1 ]
  
  # New name is the corresponding nameNoIso
  longestIsoNewName <- dt$nameNoIso[ dt$isoLengthRank==1 ]
  
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
