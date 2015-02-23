# This script generates the job scripts for the entire pipeline

pipelineOutDir <- "/mnt/NOBACKUP/mariansc/share"

species <- c( "NaSt1", "NaSt2", "MeNu", "HoVu", "StLa", "BrDi" )

readFilesTbl <- read.csv("indata/raw_read_sheet.csv")

# how to move the old files into the pipeline?


###
#
# General jobs:
# ======================
#
# Read trimming and quality control of samples

#   trimmo - trim reads
#     input: Raw reads
#     output: trimmed reads

source("processes/trimmo/createTrimmoJob.R")

createTrimmoJob(readFileFolders = readFilesTbl$PATH, outNames = readFilesTbl$sampleID,
                outDir = file.path(pipelineOutDir,"trimmo"), jobArraySize = 10, 
                adapterFile = "indata/TruSeq3-PE-2.fa")


###
#
# Genome specific jobs:
# ======================
#
#   There will be a folder for each species/de-novo assembled trancriptome. 
# Jobs that can be run on each species independantly shall reside in the 
# corresponding species folder.
#


#   trinity - De-novo transcript assembly
#     input: trimmed reads
#     output: assembled transcript sequences

#   RSEM - read counts
#     input: trimmed reads, assembled transcript sequences
#     output: read counts (genes/isoforms)

#   transDecoder - ORF finding
#     input: assembled transcript sequences
#     output: ORF sequences (pep/nucleotides)

#   longestORF - Select the longest ORF from each gene
#     input: ORF sequences (pep/nucleotides)
#     output: longest ORF (pep/nucleotides)



###
#
# Cross species jobs:
# ===================
#
# Jobs that

# orthoMCL - ortholog group finder
#   input: longetsORF sequences (pep/nucleotides)





