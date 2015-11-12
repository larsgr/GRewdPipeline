# This script generates the job scripts for the entire pipeline

readFilesTbl <- read.csv("indata/raw_read_sheet.csv", stringsAsFactors=F)

# Set umask so that new files will have the group write access
Sys.umask(mode="0002")

###
# Define directory structure
pipelineOutDir     <- "/mnt/NOBACKUP/mariansc/share"
trinityOutDir      <- file.path(pipelineOutDir,"trinity")
RSEMOutDir         <- file.path(pipelineOutDir,"RSEM") 
transdecoderOutDir <- file.path(pipelineOutDir,"transdecoder")
orthoOutDir        <- file.path(pipelineOutDir,"orthos")
refGenDir          <- file.path(pipelineOutDir,"refGenomes")

# dir.create(trinityOutDir)
# dir.create(RSEMOutDir)
# dir.create(transdecoderOutDir)
# dir.create(orthoOutDir)
# dir.create(refGenDir)


###
# source R files
source("processes/trimmo/createTrimmoJob.R")
source("processes/trinity/createTrinityJob.R")
source("processes/RSEM/createRSEMJob.R")
source("processes/transdecoder/createTransdecoderJob.R")
source("processes/filterLongestORF/createFilterORFjob.R")
source("processes/orthoMCL/createOrthoMCLjob.R")
source("processes/makeExprTables/createExprTablesJob.R")
source("processes/RJob/RJob.R")
source("processes/MAFFT/MAFFTJob.R")
source("processes/ArrayR/ArrayRJob.R")
source("processes/SLURMscript/createSLURMscript.R")


###
#
# Trimmo - Read trimming and quality control of samples
#
#     input: Raw reads
#     output: trimmed reads


createTrimmoJob(readFileFolders = readFilesTbl$PATH, outNames = readFilesTbl$sampleID,
                outDir = file.path(pipelineOutDir,"trimmo"), jobArraySize = 10, 
                adapterFile = "indata/TruSeq3-PE-2.fa")

# add the trimmed reads to the readFilesTbl
readFilesTbl$trimmedLeft <- file.path(pipelineOutDir,"trimmo",paste0(readFilesTbl$sampleID,".R1.fq"))
readFilesTbl$trimmedRight <- file.path(pipelineOutDir,"trimmo",paste0(readFilesTbl$sampleID,".R2.fq"))

###
#
#   Trinity - De-novo transcript assembly
#
#     input: trimmed reads
#     output: assembled transcript sequences

assemblies <- na.omit(unique(readFilesTbl$assembly))
trinityOutput <- list()
for(assemblyName in assemblies){
  trinityOutput[[assemblyName]] <- file.path(trinityOutDir, assemblyName, paste0(assemblyName,".fasta"))
  try( createTrinityJob(
         leftReadFiles = readFilesTbl$trimmedLeft[ readFilesTbl$assembly %in% assemblyName ],
         rightReadFiles = readFilesTbl$trimmedRight[ readFilesTbl$assembly %in% assemblyName ],
         outDir=dirname(trinityOutput[[assemblyName]]),
         trinityOutputName=basename(trinityOutput[[assemblyName]]),
         max_memory="400G", CPU=32 ) )
}


###
#
#   RSEM - read counts
#
#     input: trimmed reads, assembled transcript sequences
#     output: read counts (genes/isoforms)

# create RSEM jobs for each assembly
RSEMout <- list()
for(assemblyName in assemblies){
  idx <- readFilesTbl$assembly %in% assemblyName

  files <- file.path(RSEMOutDir,assemblyName,paste0(readFilesTbl$sampleID[idx],".genes.results"))
  names(files) <- readFilesTbl$sampleID[idx]
  RSEMout[[assemblyName]] <- files
  
  try(
    createRSEMJob( outDir = file.path(RSEMOutDir,assemblyName),
                 transcriptsFile = trinityOutput[[assemblyName]],
                 leftReadFiles = readFilesTbl$trimmedLeft[idx],
                 rightReadFiles = readFilesTbl$trimmedRight[idx],
                 outputPrefixes = readFilesTbl$sampleID[idx],
                 jobName = paste0(assemblyName,"RSEM"),
                 CPU=4, arraySize=4) )
}


###
#
#   transDecoder - ORF finding
#
#     input: assembled transcript sequences
#     output: ORF sequences (pep/nucleotides)

for( assemblyName in assemblies){
  createTransdecoderJob(outDir = file.path(transdecoderOutDir,assemblyName),
                        transcriptsFile = trinityOutput[[assemblyName]],
                        jobName = paste0(assemblyName,"TD"),
                        CPU = 1)
}

#
# Get lolium perenne transcriptome
#
# (LoPeF) Lolium Perenne low quality (genotype Falster) downloaded from:
# http://www.ebi.ac.uk/arrayexpress/files/E-MTAB-2623/E-MTAB-2623.processed.2.zip
# (LoPe) Lolium Perenne high quality downloaded from:
# http://www.ebi.ac.uk/ena/data/view/GAYX01000001-GAYX01185833

LoPe_FastaFile <- file.path(pipelineOutDir,"lolium/lolium_perenne_P226_135_16.fasta")
LoPeF_FastaFile <- file.path(pipelineOutDir,"lolium/FALSTER_transcriptome_assembly.fasta")

#
# Run transdecoder on lolium
#
createTransdecoderJob(outDir = file.path(transdecoderOutDir,"LoPeF"),
                      transcriptsFile = LoPeF_FastaFile,
                      jobName = "LoPeF_TD",
                      CPU = 1)

createTransdecoderJob(outDir = file.path(transdecoderOutDir,"LoPe"),
                      transcriptsFile = LoPe_FastaFile,
                      jobName = "LoPe_TD",
                      CPU = 1)


###
#
#   longestORF - Select the longest ORF from each gene
#
#     input: ORF sequences (pep/nucleotides)
#     output: longest ORF (pep/nucleotides)

# TODO: WARNING! This doesn't work if not transdecoder is already finished
# NOTE: Only using the MeNu1 assembly from this point
transdecoderOutPepFiles <- sapply(c("NaSt","BrDi","MeNu1","StLa","HoVu","LoPe","LoPeF"), function(x){
  dir(file.path(transdecoderOutDir,x),pattern="\\.fasta\\.transdecoder\\.pep$",full.names = T)
})


transdecoderOutCdsFiles <- sapply(c("NaSt","BrDi","MeNu1","StLa","HoVu","LoPe","LoPeF"), function(x){
  dir(file.path(transdecoderOutDir,x),pattern="\\.transdecoder\\.cds$",full.names = T)
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

refGenomes <- list(
  Bd_R = file.path(refGenDir,"brachypodium_1.2_Protein_representative.fa"),
  Hv_R = file.path(refGenDir,"barley_HighConf_genes_MIPS_23Mar12_ProteinSeq.fa")
)

refGenomesNucl <- list(
  Bd_R = file.path(refGenDir,"brachypodium_1.2_CDS.fa"),
  Hv_R = file.path(refGenDir,"barley_HighConf_genes_MIPS_23Mar12_CDSSeq.fa")
)

# download.file("ftp://ftpmips.helmholtz-muenchen.de/plants/brachypodium/v1.2/brachypodium_1.2_Protein_representative.fa",
#               destfile = refGenomes$Bd_R )
# download.file("ftp://ftpmips.helmholtz-muenchen.de/plants/barley/public_data/genes/barley_HighConf_genes_MIPS_23Mar12_ProteinSeq.fa",
#               destfile = refGenomes$Hv_R )
# download.file("ftp://ftpmips.helmholtz-muenchen.de/plants/brachypodium/v1.2/brachypodium_1.2_CDS.fa",
#               destfile = refGenomesNucl$Bd_R )
# download.file("ftp://ftpmips.helmholtz-muenchen.de/plants/barley/public_data/genes/barley_HighConf_genes_MIPS_23Mar12_CDSSeq.fa",
#               destfile = refGenomesNucl$Hv_R )

#
# Download outgroup reference species
#

outGroupGenomes <- list(
  Zm_R = file.path(refGenDir,"ZmB73_5a_WGS_translations.fasta"),
  Sb_R = file.path(refGenDir,"sorghum1.4Proteins.fa"),
  Os_R = file.path(refGenDir,"rap2Protein.fa")
)

outGroupGenomesNucl <- list(
  Zm_R = file.path(refGenDir,"ZmB73_5a_WGS_cds.fasta"),
  Sb_R = file.path(refGenDir,"sorghum1.4CDS.fa"),
  Os_R = file.path(refGenDir,"rap2BestGuessCds.fa")
)

# NOTE: Zm transcripts ID's differ for peptide ID's  (_T/_FGT instead of _P/_FGP)
# NOTE: There are multiple isoforms in these annotations..

# Dowloaded files:
# ftp://ftpmips.helmholtz-muenchen.de/plants/sorghum/sorghum1.4Proteins.fa
# ftp://ftpmips.helmholtz-muenchen.de/plants/sorghum/sorghum1.4CDS.fa
# ftp://ftp.maizesequence.org/pub/maize/release-5b/working-set/ZmB73_5a_WGS_cds.fasta.gz
# ftp://ftp.maizesequence.org/pub/maize/release-5b/working-set/ZmB73_5a_WGS_translations.fasta.gz
# ftp://ftpmips.helmholtz-muenchen.de/plants/rice/rap2BestGuessCds.fa
# ftp://ftpmips.helmholtz-muenchen.de/plants/rice/rap2Protein.fa



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


createOrthoMCLjob( outDir = file.path(orthoOutDir,"orthoMCL"),
                   proteomeFiles = c(filterORFout,
                                     unlist(refGenomes),
                                     unlist(outGroupGenomes)),
                   taxon_codes = c(names(transdecoderOutPepFiles),
                                   names(refGenomes),
                                   names(outGroupGenomes)),
                   blastCPU=1, blastArraySize=150)

orthoMCLout <- file.path(orthoOutDir,"orthoMCL","groups.txt")

# makeExprTables - Convert RSEM output files to more handy expression tables
#   input: groups file from orthoMCL
#          Expression counts files from RSEM
#   output: Several tables

# NOTE: MeNu2 is not included anymore
RSEMoutSubset <- RSEMout
RSEMoutSubset$MeNu2 <- NULL
RSEMoutSubset$NaSt <- RSEMout$NaSt[grepl("NaSt1",names(RSEMout$NaSt))]

createExprTablesJob(outDir = file.path(orthoOutDir,"exprTbls"),
                    orthoGrpFile=orthoMCLout,
                    RSEMout=RSEMoutSubset)


###
#
# Generate phylegentic trees for each ortho grp
# ==============================================
#
# 1. Generate fasta files for each orthogroup (peptides)
# 2. Generate nucleotide fasta files for each ortholog group
# 3. Align each group ufing MAFFT (peptides)
# 4. Convert protein alignments to nucleotide alignments with Pal2Nal
# 5. Generate trees for the nucleotide alignments


#
# 1. generate fasta files for each orthogroup (peptides)
#


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


generateScript(grpFastasJob)

#
# 2. Generate nucleotide fasta files for each ortholog group
#

# Get corresponding nucleotide sequences by looking up in the tables
# generated when selecting the longest ORF's, or use pepID2nucID function (for ref genomes)

RJob( outDir = file.path(orthoOutDir,"grpCDSFastas"),jobName = "grpCDSFastas",
      data = list( orthoGrpFile = orthoMCLout, 
                   cdsFastas = c(refGenomesNucl,outGroupGenomesNucl,transdecoderOutCdsFiles),
                   filterORFtbl.out = 
                     setNames(file.path(orthoOutDir,"longestORFs",
                                        paste0(names(transdecoderOutPepFiles),".tbl")),
                              names(transdecoderOutPepFiles)),
                   pepID2nucID = list(
                     Zm_R = function(pepID){
                       sub( "_FGP","_FGT",
                            sub("_P","_T",pepID))
                     }) ),
      FUN = function(){
        source("/mnt/users/lagr/GRewd/pipeline/R/orthoGrpTools.R")
        # load ortholog groups
        grps <- loadOrthoGrpsArray(orthoGrpFile = data$orthoGrpFile)
        
        # load sequence ID conversion table
        lapply(data$filterORFtbl.out, function(tblFile) {
          tmp <- readr::read_tsv(tblFile,col_names = c("seqID","cdsID"))
          setNames(tmp$cdsID,tmp$seqID)
        })  -> cdsID
        
        # load all input sequences
        lapply(data$cdsFastas, function( faFile ) {
          seqinr::read.fasta(faFile)
        }) -> seqs
        
        # for each ortho group:
        for( i in 1:nrow(grps)){
          
          fullSeqIDs <- grpToChar(grps[i,]) # names of the sequences
          
          if(length(fullSeqIDs)>4){ # only bother with groups of more than 4
            
            outSeqs <- list() # list to be filled sequences
            
            # for each sequence in group:
            for(fullSeqID in fullSeqIDs){
              spc <- sub("\\|.*","",fullSeqID)
              seqID <- sub(".*\\|","",fullSeqID)
              if(spc %in% names(cdsID)){
                outSeqs <- c(outSeqs,seqs[[spc]][ cdsID[[spc]][seqID] ])
              } else if(spc %in% names(data$pepID2nucID)){
                outSeqs <- c(outSeqs,seqs[[spc]][ data$pepID2nucID[[spc]](seqID) ])
              } else {
                outSeqs <- c(outSeqs,seqs[[spc]][ seqID ])
              }
            }
            
            # write sequences to file:
            seqinr::write.fasta(outSeqs,names=fullSeqIDs,
                                file.out = paste0(rownames(grps)[i],".cds"))
          }
        }
      }) -> grpCDSFastasJob

generateScript(grpCDSFastasJob)

#
# 3. align each group ufing MAFFT (peptides)
#

#
MAFFTJob(outDir = file.path(orthoOutDir,"grpAligned"),
         arraySize = 100,
         inFastaDir = grpFastasJob$outDir 
         ) -> myMAFFTJob

generateScript(myMAFFTJob)



#
# 4. Convert protein alignments to nucleotide alignments with Pal2Nal
# 




# for each alignment (with more than four sequences)
arraySize = 100
ArrayRJob(x = 1:arraySize, 
          outDir = file.path(orthoOutDir,"pal2nal"),
          jobName = "pal2nal",
          commonData=list( grpCDSpath = grpCDSFastasJob$outDir,
                           grpAlignedPepPath = myMAFFTJob$outDir,
                           orthoGrpFile = orthoMCLout,
                           arraySize = arraySize,
                           pipelineSrcPath = getwd()),
          FUN=function(x){
            source(file.path(commonData$pipelineSrcPath,"R/orthoGrpTools.R"))
            
            
            grpTbl <- loadOrthoGrpsTable(orthoGrpFile = commonData$orthoGrpFile)
            grpSizes <- table(grpTbl$grpID)
            # only use the groups with more than four sequences
            alnBigFaFiles <- file.path( commonData$grpAlignedPepPath,
                                        paste0(names(which(grpSizes>4)),".aln"))
            
            for( i in seq(x,length(alnBigFaFiles),by = commonData$arraySize)){
              alignedPepFile <- alnBigFaFiles[i]
              cdsFile <- file.path(commonData$grpCDSpath, 
                                   sub("aln$","cds",basename(alignedPepFile)))
              alignedCdsFile <- sub("aln$","cds.aln",basename(alignedPepFile))
              system(paste(sep="\n",
                           "module load pal2nal/14.0",
                           "module load mafft/7.130",
                           "",
                           paste("pal2nal.pl",
                                 alignedPepFile,
                                 cdsFile,
                                 "-output fasta",
                                 ">", alignedCdsFile) ))
            }
          }) -> pal2nalJob

generateScript(pal2nalJob)


#
# 5. Generate trees for the nucleotide alignments
#

arraySize = 1000
ArrayRJob(x = 1:arraySize, 
          outDir = file.path(orthoOutDir,"treesNuc"),
          jobName = "nucTree", 
          commonData=list( grpAlignedCdsPath = pal2nalJob$outDir,
                           arraySize = arraySize,
                           pipelineSrcPath = getwd()),
          FUN=function(x){
            source(file.path(commonData$pipelineSrcPath,"processes/phangorn/makeTree.R"))
            
            grpAlignedCdsFiles <- dir(commonData$grpAlignedCdsPath, full.names = T,
                                      pattern="\\.aln$" )
            
            grpAlignedCdsFiles <- rev(grpAlignedCdsFiles) # smallest grps first
            
            for( i in seq(x,length(grpAlignedCdsFiles),by = commonData$arraySize,)){
              alignedFastaFile <- grpAlignedCdsFiles[i]
              
              makeTree(alignedFastaFile = alignedFastaFile, outDir = ".", type="DNA")
            }
          }) -> makeNucTreeJob

generateScript(makeNucTreeJob)




##
#
# Split complex trees with clanFinder:
#
# Output:
#   splitGroups.txt - New (sub)groups resulting from split trees
#   goodTrees.rds - trees that pass minimum requirement (including split trees)
#   goodGroups.txt - corresponding grps of the goodTrees
#   goodTreeStats.rds - table with stats about each tree

RJob(jobName = "splitGroups", outDir = file.path(orthoOutDir,"splitGroups"),
     data = list( treePath = makeNucTreeJob$outDir,
                  outGrpFile = "splitGroups.txt"),
     FUN=function(){
       source("/mnt/users/lagr/GRewd/pipeline/processes/splitTreesToGrps/splitTreesToGrps.R")
       with(data,
         splitTreesToGrps(treePath)
       )
     }) -> splitGroupsJob

generateScript(splitGroupsJob)



####
# Run DESeq
#
# Input: exprTbls
# Output:
#   DE.RData - DE object containing vst, resRamp and resPeak for each species.

dir.create(file.path(orthoOutDir,"DESeq"))
createSLURMscript( jobName = "runDESeq", workdir = file.path(orthoOutDir,"DESeq"),
                    script = c("module load R/3.1.0",
                               "Rscript /mnt/users/lagr/GRewd/pipeline/processes/DESeq/runDESeq.R")
                  )

#########
#
# Run codeml
#
arraySize = 100
ArrayRJob( x = 1:arraySize, outDir = file.path(orthoOutDir,"PAML"),
           jobName = "codeml", arraySize = arraySize, 
           commonData = list(
             goodTreesFile = "/mnt/NOBACKUP/mariansc/share/orthos/splitGroups/goodTrees.rds",
             goodTreeStatFile = "/mnt/NOBACKUP/mariansc/share/orthos/splitGroups/goodTreeStats.rds",
             codonAlnPath = "/mnt/NOBACKUP/mariansc/share/orthos/pal2nal",
             arraySize = arraySize),
           FUN = function(x){
             source("/mnt/users/lagr/GRewd/pipeline/processes/runPAML/runPAML.R")
             mainLoopPAML(x, 
                          goodTreesFile = commonData$goodTreesFile, 
                          goodTreeStatFile = commonData$goodTreeStatFile,
                          codonAlnPath = commonData$codonAlnPath,
                          arraySize = commonData$arraySize )
           }) -> PAMLjob

generateScript(PAMLjob)




