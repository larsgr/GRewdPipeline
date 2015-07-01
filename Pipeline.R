# This script generates the job scripts for the entire pipeline

pipelineOutDir <- "/mnt/NOBACKUP/mariansc/share"

readFilesTbl <- read.csv("indata/raw_read_sheet.csv", stringsAsFactors=F)


# Set umask so that new files will have the group write access
Sys.umask(mode="0002")


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

# BrDi 7 samples (4 correct and 3 wrong)
asmSamples$BrDi <- grepl("BrDi",readFilesTbl$SPECIES) & readFilesTbl$chemistry == "correct"
asmSamples$wc_BrDi <- grepl("BrDi",readFilesTbl$SPECIES) & readFilesTbl$chemistry == "wrong"

# HoVu 11 samples (8 correct and 3 wrong)
asmSamples$HoVu <- grepl("HoVu",readFilesTbl$SPECIES) & readFilesTbl$chemistry == "correct"
asmSamples$wc_HoVu <- grepl("HoVu",readFilesTbl$SPECIES) & readFilesTbl$chemistry == "wrong"
# MeNu1 37 samples (21 correct and 16 wrong)
asmSamples$MeNu1 <- grepl("MeNu1",readFilesTbl$SPECIES) & readFilesTbl$chemistry == "correct"
asmSamples$wc_MeNu1 <- grepl("MeNu1",readFilesTbl$SPECIES) & readFilesTbl$chemistry == "wrong"

# MeNu2 28 samples (16 correct and 12 wrong)
asmSamples$MeNu2 <- grepl("MeNu2",readFilesTbl$SPECIES) & readFilesTbl$chemistry == "correct"
asmSamples$wc_MeNu2 <- grepl("MeNu2",readFilesTbl$SPECIES) & readFilesTbl$chemistry == "wrong"

# (NaSt1+2+3) 27 samples (18 correct and 9 wrong)
asmSamples$NaSt <- grepl("NaSt",readFilesTbl$SPECIES) & readFilesTbl$chemistry == "correct"
asmSamples$wc_NaSt <- grepl("NaSt",readFilesTbl$SPECIES) & readFilesTbl$chemistry == "wrong"

# StLa 10 samples (7 correct and 3 wrong)
asmSamples$StLa <- grepl("StLa",readFilesTbl$SPECIES) & readFilesTbl$chemistry == "correct"
asmSamples$wc_StLa <- grepl("StLa",readFilesTbl$SPECIES) & readFilesTbl$chemistry == "wrong"


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
RSEMsamplesIdx$NaSt <- grepl("NaSt",readFilesTbl$SPECIES)
RSEMsamplesIdx$HoVu <- grepl("HoVu",readFilesTbl$SPECIES)
RSEMsamplesIdx$MeNu1 <- grepl("MeNu1",readFilesTbl$SPECIES) & readFilesTbl$chemistry == "correct"
RSEMsamplesIdx$MeNu2 <- grepl("MeNu2",readFilesTbl$SPECIES) & readFilesTbl$chemistry == "correct"
RSEMsamplesIdx$BrDi <- grepl("BrDi",readFilesTbl$SPECIES)
RSEMsamplesIdx$StLa <- grepl("StLa",readFilesTbl$SPECIES)
RSEMsamplesIdx$wc_MeNu2 <- grepl("MeNu2",readFilesTbl$SPECIES) & readFilesTbl$chemistry == "wrong"


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

RSEMout <- list()
for(assemblyName in names(RSEMsamplesIdx)){
  idx <- RSEMsamplesIdx[[assemblyName]]
  files <- file.path(RSEMOutDir,assemblyName,paste0(readFilesTbl$sampleID[idx],".genes.results"))
  names(files) <- readFilesTbl$sampleID[idx]
  RSEMout[[assemblyName]] <- files
}

#   transDecoder - ORF finding
#     input: assembled transcript sequences
#     output: ORF sequences (pep/nucleotides)


transdecoderOutDir <- file.path(pipelineOutDir,"transdecoder")
dir.create(transdecoderOutDir)

source("processes/transdecoder/createTransdecoderJob.R")

for( assemblyName in c("NaSt","BrDi","MeNu1","MeNu2","StLa","HoVu")){
  createTransdecoderJob(outDir = file.path(transdecoderOutDir,assemblyName),
                        transcriptsFile = trinityOutput[[assemblyName]],
                        jobName = paste0(assemblyName,"TD"),
                        CPU = 1)
}

#   longestORF - Select the longest ORF from each gene
#     input: ORF sequences (pep/nucleotides)
#     output: longest ORF (pep/nucleotides)

# create one directory to contain cross species ortholog related jobs

orthoOutDir <- file.path(pipelineOutDir,"orthos")
dir.create(orthoOutDir)

source("processes/filterLongestORF/createFilterORFjob.R")

transdecoderOutPepFiles <- sapply(c("NaSt","BrDi","MeNu1","MeNu2","StLa","HoVu"), function(x){
  file.path(transdecoderOutDir,x,paste0(x,".fasta.transdecoder.pep"))
})
                  
createFilterORFjob(outDir = file.path(orthoOutDir,"longestORFs"),
                   ORFfiles = transdecoderOutPepFiles,
                   outPrefix = names(transdecoderOutPepFiles))

filterORFout <- file.path(orthoOutDir,"longestORFs",paste0(names(transdecoderOutPepFiles),".longest.fasta"))


###############
#
# Download reference protein sequences for Barley and Brachypodium
#
# Note that the proteomes have only one representative sequence per gene.

refGenDir <- file.path(pipelineOutDir,"refGenomes")
refGenomes <- list(
  Bd_R = file.path(refGenDir,"brachypodium_1.2_Protein_representative.fa"),
  Hv_R = file.path(refGenDir,"barley_HighConf_genes_MIPS_23Mar12_ProteinSeq.fa")
)

refGenomesNucl <- list(
  Bd_R = file.path(refGenDir,"brachypodium_1.2_CDS.fa"),
  Hv_R = file.path(refGenDir,"barley_HighConf_genes_MIPS_23Mar12_CDSSeq.fa")
)

# dir.create(refGenDir)
# download.file("ftp://ftpmips.helmholtz-muenchen.de/plants/brachypodium/v1.2/brachypodium_1.2_Protein_representative.fa",
#               destfile = refGenomes$Bd_R )
# download.file("ftp://ftpmips.helmholtz-muenchen.de/plants/barley/public_data/genes/barley_HighConf_genes_MIPS_23Mar12_ProteinSeq.fa",
#               destfile = refGenomes$Hv_R )
# download.file("ftp://ftpmips.helmholtz-muenchen.de/plants/brachypodium/v1.2/brachypodium_1.2_CDS.fa",
#               destfile = refGenomesNucl$Bd_R )
# download.file("ftp://ftpmips.helmholtz-muenchen.de/plants/barley/public_data/genes/barley_HighConf_genes_MIPS_23Mar12_CDSSeq.fa",
#               destfile = refGenomesNucl$Hv_R )


###
#
# Cross species jobs:
# ===================
#


# orthoMCL - ortholog group finder
#   input: longetsORF sequences (pep/nucleotides)
#  output: groups.txt (ortholog groups)
#          allProteomes.fasta (combined sequences from all species)
#          allProteomes.db* (blast database)
#          all_vs_all.out (all vs all blast result)

source("processes/orthoMCL/createOrthoMCLjob.R")

createOrthoMCLjob( outDir = file.path(orthoOutDir,"orthoMCL"),
                   proteomeFiles = c(filterORFout,unlist(refGenomes)),
                   taxon_codes = c(names(transdecoderOutPepFiles),names(refGenomes)),
                   blastCPU=66)

orthoMCLout <- file.path(orthoOutDir,"orthoMCL","groups.txt")
## TODO: blast only used 7% of CPU!! Should divide the blast job into parts and execute distributed


# makeExprTables - Convert RSEM output files to more handy expression tables
#   input: groups file from orthoMCL
#          Expression counts files from RSEM
#   output: Several tables

source("processes/makeExprTables/createExprTablesJob.R")

createExprTablesJob(outDir = file.path(orthoOutDir,"exprTbls"),
                    orthoGrpFile=orthoMCLout,
                    RSEMout=RSEMout)


###
#
# Generate phylegentic trees for each ortho grp
# ==============================================
#

# 1. generate fasta files for each orthogroup

source("processes/RJob/RJob.R")

RJob( outDir = file.path(orthoOutDir,"grpFastas"),
      data = list( orthoGrpFile = orthoMCLout, 
                   inFasta = file.path(orthoOutDir,"orthoMCL","allProteomes.fasta")), 
      FUN= function(){
        source("/mnt/users/lagr/GRewd/pipeline/R/orthoGrpTools.R")
        grps <- loadOrthoGrpsArray(orthoGrpFile = data$orthoGrpFile)
        seqs <- seqinr::read.fasta(file = data$inFasta)
        for( i in 1:nrow(grps)){
          seqIDs <- grpToChar(grps[i,])
          seqinr::write.fasta(seqs[seqIDs],names=seqIDs,
                              file.out = paste0(rownames(grps)[i],".fasta"))
        }
     }) -> grpFastasJob


# generateScript(grpFastasJob)

# 2. align each group ufing MAFFT

source("processes/MAFFT/MAFFTJob.R")

#
MAFFTJob(outDir = file.path(orthoOutDir,"grpAligned"),
         arraySize = 100,
         inFastaDir = grpFastasJob$outDir 
         ) -> myMAFFTJob

# generateScript(myMAFFTJob)


# 3. Make trees using phangorn
#
# TODO: Figuring out which groups to use should be done as part of the job

source("processes/ArrayR/ArrayRJob.R")

source("/mnt/users/lagr/GRewd/pipeline/R/orthoGrpTools.R")
grpTbl <- loadOrthoGrpsTable(orthoGrpFile = orthoMCLout)
grpSizes <- table(grpTbl$grpID)
# only use the groups with more than four sequences
alnBigFaFiles <- file.path(myMAFFTJob$outDir,paste0(names(which(grpSizes>4)),".aln"))

ArrayRJob(x = rev(alnBigFaFiles), outDir = file.path(orthoOutDir,"trees"),
          jobName = "makeTree", arraySize = 100, 
          FUN=function(alignedFastaFile){
            source("/mnt/users/lagr/GRewd/pipeline/processes/phangorn/makeTree.R")
            makeTree(alignedFastaFile = alignedFastaFile, outDir = ".")
          }) -> makeTreeJob

generateScript(makeTreeJob)
