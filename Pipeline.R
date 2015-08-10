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


                                                       

source("processes/trinity/createTrinityJob.R")

assemblies <- na.omit(unique(readFilesTbl$assembly))

trinityOutput <- list()
for(assemblyName in assemblies){
  trinityOutput[[assemblyName]] <- file.path(pipelineOutDir,"trinity", assemblyName, paste0(assemblyName,".fasta"))
  try( createTrinityJob(
         leftReadFiles = readFilesTbl$trimmedLeft[ readFilesTbl$assembly %in% assemblyName ],
         rightReadFiles = readFilesTbl$trimmedRight[ readFilesTbl$assembly %in% assemblyName ],
         outDir=dirname(trinityOutput[[assemblyName]]),
         trinityOutputName=basename(trinityOutput[[assemblyName]]),
         max_memory="400G", CPU=32 ) )
}


#   RSEM - read counts
#     input: trimmed reads, assembled transcript sequences
#     output: read counts (genes/isoforms)

source("processes/RSEM/createRSEMJob.R")


# Put all read count jobs in sub-folders of RSEM
RSEMOutDir <- file.path(pipelineOutDir,"RSEM")
dir.create(RSEMOutDir)


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


#   transDecoder - ORF finding
#     input: assembled transcript sequences
#     output: ORF sequences (pep/nucleotides)


transdecoderOutDir <- file.path(pipelineOutDir,"transdecoder")
dir.create(transdecoderOutDir)

source("processes/transdecoder/createTransdecoderJob.R")

for( assemblyName in assemblies){
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

transdecoderOutCdsFiles <- sapply(c("NaSt","BrDi","MeNu1","MeNu2","StLa","HoVu"), function(x){
  file.path(transdecoderOutDir,x,paste0(x,".fasta.transdecoder.cds"))
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

source("processes/orthoMCL/createOrthoMCLjob.R")

createOrthoMCLjob( outDir = file.path(orthoOutDir,"orthoMCL"),
                   proteomeFiles = c(filterORFout,
                                     unlist(refGenomes),
                                     unlist(outGroupGenomes)),
                   taxon_codes = c(names(transdecoderOutPepFiles),
                                   names(refGenomes),
                                   names(outGroupGenomes)),
                   blastCPU=1, blastArraySize=100)

orthoMCLout <- file.path(orthoOutDir,"orthoMCL","groups.txt")

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
# 1. Generate fasta files for each orthogroup (peptides)
# 2. Generate nucleotide fasta files for each ortholog group
# 3. Align each group ufing MAFFT (peptides)
# 4. Convert protein alignments to nucleotide alignments with Pal2Nal
# 5. Generate trees for the nucleotide alignments


#
# 1. generate fasta files for each orthogroup (peptides)
#

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

source("processes/MAFFT/MAFFTJob.R")

#
MAFFTJob(outDir = file.path(orthoOutDir,"grpAligned"),
         arraySize = 100,
         inFastaDir = grpFastasJob$outDir 
         ) -> myMAFFTJob

generateScript(myMAFFTJob)



#
# 4. Convert protein alignments to nucleotide alignments with Pal2Nal
# 

source("processes/ArrayR/ArrayRJob.R")


# TODO: Figuring out which groups to use should be done as part of the job
source("R/orthoGrpTools.R")
grpTbl <- loadOrthoGrpsTable(orthoGrpFile = orthoMCLout)
grpSizes <- table(grpTbl$grpID)
# only use the groups with more than four sequences
alnBigFaFiles <- file.path(myMAFFTJob$outDir,paste0(names(which(grpSizes>4)),".aln"))

# for each alignment (with more than four sequences)
ArrayRJob(x = rev(alnBigFaFiles), outDir = file.path(orthoOutDir,"pal2nal"),
          jobName = "pal2nal", arraySize = 10,
          commonData=list( grpCDSpath = grpCDSFastasJob$outDir),
          FUN=function(alignedPepFile){
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
          }) -> pal2nalJob

generateScript(pal2nalJob)

alnBigCdsFiles <- file.path(pal2nalJob$outDir,paste0(names(which(grpSizes>4)),".cds.aln"))

#
# 5. Generate trees for the nucleotide alignments
#

ArrayRJob(x = rev(alnBigCdsFiles), outDir = file.path(orthoOutDir,"treesNuc"),
          jobName = "nucTree", arraySize = 100, 
          FUN=function(alignedFastaFile){
            source("/mnt/users/lagr/GRewd/pipeline/processes/phangorn/makeTree.R")
            makeTree(alignedFastaFile = alignedFastaFile, outDir = ".", type="DNA")
          }) -> makeNucTreeJob

generateScript(makeNucTreeJob)




##
#
# Split complex trees with clanFinder:
#

splitOrthosDir <- file.path(orthoOutDir,"splitOrthos")
dir.create(splitOrthosDir)

RJob(jobName = "clanFinder", outDir = file.path(splitOrthosDir,"clanFinder"),
     data = list( treePath = makeNucTreeJob$outDir,
                  outGrpFile = "splitGroups.txt"),
     FUN=function(){
       source("/mnt/users/lagr/GRewd/pipeline/processes/splitTreesToGrps/splitTreesToGrps.R")
       with(data,
         splitTreesToGrps(treePath, outGrpFile, outSpcs = c("Os_R", "Sb_R", "Zm_R"))
       )
     }) -> clanFinderJob

generateScript(clanFinderJob)

splitGrpFile <- file.path(clanFinderJob$outDir,clanFinderJob$params$data$outGrpFile)

###
#
# Generate phylegentic trees for each split ortho grp
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


RJob( outDir = file.path(splitOrthosDir,"grpFastas"), jobName = "splitGrpFastas",
      data = list( orthoGrpFile = splitGrpFile, 
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
      }) -> splitGrpFastasJob


generateScript(splitGrpFastasJob)

#
# 2. Generate nucleotide fasta files for each ortholog group
#

# Get corresponding nucleotide sequences by looking up in the tables
# generated when selecting the longest ORF's, or use pepID2nucID function (for ref genomes)

RJob( outDir = file.path(splitOrthosDir,"grpCDSFastas"),jobName = "splitGrpCDSFastas",
      data = list( orthoGrpFile = splitGrpFile, 
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
      }) -> splitGrpCDSFastasJob

generateScript(splitGrpCDSFastasJob)

#
# 3. align each group ufing MAFFT (peptides)
#

MAFFTJob(outDir = file.path(splitOrthosDir,"grpAligned"), jobName = "splitMAFFT",
         arraySize = 10,
         inFastaDir = splitGrpFastasJob$outDir
) -> splitMAFFTJob

generateScript(splitMAFFTJob)



#
# 4. Convert protein alignments to nucleotide alignments with Pal2Nal
# 


N <- 10 # arraySize
ArrayRJob(x = 1:N, outDir = file.path(splitOrthosDir,"pal2nal"),
          jobName = "splitPal2nal",
          commonData=list( grpCDSpath = splitGrpCDSFastasJob$outDir,
                           alignedPepPath = splitMAFFTJob$outDir,
                           N = N),
          FUN=function(x){
            # list aligned pep files
            alignedPepFiles <- dir(commonData$alignedPepPath,pattern="aln$",full.names = T)
            
            # run mafft on every N'th file, starting with file x
            for(i in seq(x,length(alignedPepFiles),by=commonData$N)){
              alignedPepFile <- alignedPepFiles[i]
              
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
          }) -> splitPal2nalJob

generateScript(splitPal2nalJob)


#
# 5. Generate trees for the nucleotide alignments
#

N <- 50 # arraySize
ArrayRJob(x = 1:N, outDir = file.path(splitOrthosDir,"treesNuc"),
          jobName = "splitNucTree", 
          commonData=list( alignedCdsPath = splitPal2nalJob$outDir,
                           N = N),
          FUN=function(x){
            source("/mnt/users/lagr/GRewd/pipeline/processes/phangorn/makeTree.R")

            # list aligned cds files
            alignedCdsFiles <- dir(commonData$alignedCdsPath,pattern="aln$",full.names = T)
            
            # run mafft on every N'th file, starting with file x
            failedTrees <- list()
            for(i in seq(x,length(alignedCdsFiles),by=commonData$N)){
              
              cat("Generating tree for", alignedCdsFiles[i],"\n")
              
              tryCatch(
                makeTree(alignedFastaFile = alignedCdsFiles[i], outDir = ".", type="DNA", bootstrap=0),
                error = function(e) {
                  failedTrees <- c(failedTrees,alignedCdsFiles[i])
                  cat("Error occured when generating tree for", alignedCdsFiles[i],":",e$message,"\n")
                })
              # throw error if any of the trees failed
              if(length(failedTrees)>0){
                stop( paste("Number of failed trees:",length(failedTrees)))
              }            }

          }) -> makeSplitNucTreeJob

generateScript(makeSplitNucTreeJob)

