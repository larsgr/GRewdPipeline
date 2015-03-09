# This script generates the job scripts for the entire pipeline

pipelineOutDir <- "/mnt/NOBACKUP/mariansc/share"

species <- c( "NaSt1", "NaSt2", "MeNu", "HoVu", "StLa", "BrDi" )

readFilesTbl <- read.csv("indata/raw_read_sheet.csv", stringsAsFactors=F)

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

# Which samples to include in each assembly:

asmSamples <- list()

# BrDi 3 samples (3 correct and 3 wrong)
asmSamples$BrDi <- grepl("BrDi",readFilesTbl$SPECIES) & readFilesTbl$chemistry == "correct"
asmSamples$wc_BrDi <- grepl("BrDi",readFilesTbl$SPECIES) & readFilesTbl$chemistry == "wrong"

# HoVu 7 samples
asmSamples$HoVu <- grepl("HoVu",readFilesTbl$SPECIES)

# MeNu1 16 samples (16 correct and 16 wrong)
asmSamples$MeNu1 <- grepl("MeNu1",readFilesTbl$SPECIES) & readFilesTbl$chemistry == "correct"
asmSamples$wc_MeNu1 <- grepl("MeNu1",readFilesTbl$SPECIES) & readFilesTbl$chemistry == "wrong"

# MeNu2 12 samples
asmSamples$MeNu2 <- grepl("MeNu2",readFilesTbl$SPECIES)

# (NaSt1+2+3) 17 samples (just correct)
asmSamples$NaSt <- grepl("NaSt",readFilesTbl$SPECIES) & readFilesTbl$chemistry == "correct"

# StLa 7 samples (just correct)
asmSamples$StLa <- grepl("StLa",readFilesTbl$SPECIES) & readFilesTbl$chemistry == "correct"


# define output files
lapply( setNames(names(asmSamples), names(asmSamples)),function(assemblyName){
  file.path(pipelineOutDir,
            paste0("trinity_",assemblyName),
            paste0(assemblyName,".fasta"))
}) -> trinityOutput
                                                         

source("processes/trinity/createTrinityJob.R")


for(assemblyName in names(asmSamples)){
  createTrinityJob(leftReadFiles = readFilesTbl$trimmedLeft[asmSamples[[assemblyName]]],
                   rightReadFiles = readFilesTbl$trimmedRight[asmSamples[[assemblyName]]],
                   outDir=dirname(trinityOutput[[assemblyName]]),
                   trinityOutputName=basename(trinityOutput[[assemblyName]]),
                   max_memory="400G", CPU=32)
}


#   RSEM - read counts
#     input: trimmed reads, assembled transcript sequences
#     output: read counts (genes/isoforms)

source("processes/RSEM/createRSEMJob.R")


# Put all read count jobs in sub-folders of RSEM
RSEMOutDir <- file.path(pipelineOutDir,"RSEM")
dir.create(RSEMOutDir)


# Define which samples to map to which assembly
RSEMsamplesIdx <- list()
RSEMsamplesIdx$HoVu <- grepl("HoVu",readFilesTbl$SPECIES)
RSEMsamplesIdx$MeNu1 <- grepl("MeNu1",readFilesTbl$SPECIES)
RSEMsamplesIdx$MeNu2 <- grepl("MeNu2",readFilesTbl$SPECIES)
RSEMsamplesIdx$BrDi <- grepl("BrDi",readFilesTbl$SPECIES)
RSEMsamplesIdx$StLa <- grepl("StLa",readFilesTbl$SPECIES)

# create RSEM job for HoVu samples
for(assemblyName in names(RSEMsamplesIdx)){
  idx <- RSEMsamplesIdx[[assemblyName]]
  createRSEMJob( outDir = file.path(RSEMOutDir,assemblyName),
                 transcriptsFile = trinityOutput[[assemblyName]],
                 leftReadFiles = readFilesTbl$trimmedLeft[idx],
                 rightReadFiles = readFilesTbl$trimmedRight[idx],
                 outputPrefixes = readFilesTbl$sampleID[idx],
                 jobName = paste0(assemblyName,"RSEM"),
                 CPU=4, arraySize=4)  
}


#   transDecoder - ORF finding
#     input: assembled transcript sequences
#     output: ORF sequences (pep/nucleotides)


transdecoderOutDir <- file.path(pipelineOutDir,"transdecoder")
dir.create(transdecoderOutDir)

source("processes/transdecoder/createTransdecoderJob.R")

assemblyName="HoVu"
createTransdecoderJob(outDir = file.path(transdecoderOutDir,assemblyName),
                      transcriptsFile = trinityOutput[[assemblyName]],
                      jobName = paste0(assemblyName,"TD"),
                      CPU = 60)

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





