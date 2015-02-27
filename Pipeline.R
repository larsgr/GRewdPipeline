# This script generates the job scripts for the entire pipeline

pipelineOutDir <- "/mnt/NOBACKUP/mariansc/share"

species <- c( "NaSt1", "NaSt2", "MeNu", "HoVu", "StLa", "BrDi" )

readFilesTbl <- read.csv("indata/raw_read_sheet.csv")

# Set umask so that new files will have the group write access
Sys.umask(mode="0002")

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

# add the trimmed reads to the readFilesTbl
readFilesTbl$trimmedLeft <- file.path(pipelineOutDir,"trimmo",paste0(readFilesTbl$sampleID,".R1.fq"))
readFilesTbl$trimmedRight <- file.path(pipelineOutDir,"trimmo",paste0(readFilesTbl$sampleID,".R2.fq"))

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


# assemblies:
# BrDi 3 samples
# HoVu 7 samples
# MeNu1 16 samples
# MeNu2 12 samples
# (NaSt1+2+3) 17 samples
# StLa 7 samples

source("processes/trinity/createTrinityJob.R")

createTrinityJob(leftReadFiles = readFilesTbl$trimmedLeft[grepl("BrDi",readFilesTbl$SPECIES)],
                 rightReadFiles = readFilesTbl$trimmedRight[grepl("BrDi",readFilesTbl$SPECIES)],
                 outDir=file.path(pipelineOutDir,"trinity_BrDi"), max_memory="400G", CPU=32)

createTrinityJob(leftReadFiles = readFilesTbl$trimmedLeft[grepl("MeNu1",readFilesTbl$SPECIES)],
                 rightReadFiles = readFilesTbl$trimmedRight[grepl("MeNu1",readFilesTbl$SPECIES)],
                 outDir=file.path(pipelineOutDir,"trinity_MeNu1"), max_memory="400G", CPU=32)

createTrinityJob(leftReadFiles = readFilesTbl$trimmedLeft[grepl("MeNu2",readFilesTbl$SPECIES)],
                 rightReadFiles = readFilesTbl$trimmedRight[grepl("MeNu2",readFilesTbl$SPECIES)],
                 outDir=file.path(pipelineOutDir,"trinity_MeNu2"), max_memory="400G", CPU=32)

createTrinityJob(leftReadFiles = readFilesTbl$trimmedLeft[grepl("NaSt",readFilesTbl$SPECIES)],
                 rightReadFiles = readFilesTbl$trimmedRight[grepl("NaSt",readFilesTbl$SPECIES)],
                 outDir=file.path(pipelineOutDir,"trinity_NaSt"), max_memory="400G", CPU=32)

createTrinityJob(leftReadFiles = readFilesTbl$trimmedLeft[grepl("StLa",readFilesTbl$SPECIES)],
                 rightReadFiles = readFilesTbl$trimmedRight[grepl("StLa",readFilesTbl$SPECIES)],
                 outDir=file.path(pipelineOutDir,"trinity_StLa"), max_memory="400G", CPU=32)

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





