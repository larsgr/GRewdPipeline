#
# Prepare files and sample attributes for submission to ArrayExpress
#

outDir <- "~/GRewd/submitArrayExpress"

readFilesTbl <- read.csv("indata/raw_read_sheet.csv", stringsAsFactors=F)


# remove files that are not supposed to be included
readFilesTbl <- readFilesTbl[!is.na(readFilesTbl$assembly),]
readFilesTbl <- readFilesTbl[!(readFilesTbl$SPECIES %in% c("NaSt2","NaSt3","MeNu2")),]

# rename samples
renameSamples <- function(sampleNames){
  for( substitusion in list( c("NaSt1","NaSt"),
                             c("T-1","W0"),
                             c("T0","D0"),
                             c("T1","D1"),
                             c("T3","W4"),
                             c("T4","W9"))){
    sampleNames <- sub(substitusion[1],substitusion[2],sampleNames)
  }
  return(sampleNames)
}

readFilesTbl$oldSampleID <- readFilesTbl$sampleID # save original sample names for later

readFilesTbl$sampleID <- renameSamples(readFilesTbl$sampleID)

# reorder to have same order as in the supplementary data
supSampleIDs <- colnames(read.table("SupplementaryData/HCOGreadCounts.tsv",header=T,nrows = 1,row.names = 1))

readFilesTbl <- readFilesTbl[match(supSampleIDs,readFilesTbl$sampleID), ]

# check that the sample names correspond to the ones in the supplementary data
identical(readFilesTbl$sampleID,supSampleIDs)


# get the path to corresponding .tar files
readFilesTbl$tarFile <- paste0(dirname(readFilesTbl$PATH),".tar")

# get fastq files in tar files
readFilesTbl$fastqR1 <- ""
readFilesTbl$fastqR2 <- ""
for(tarFile in unique(readFilesTbl$tarFile)){
  filesInTar <- system(paste("tar -tf",tarFile, "| grep fastq.gz"),intern = T)
  filesInTarR1 <- filesInTar[grepl("R1",basename(filesInTar))]
  filesInTarR2 <- filesInTar[grepl("R2",basename(filesInTar))]
  
  idx <- match( basename(readFilesTbl$PATH[readFilesTbl$tarFile==tarFile]),
                basename(dirname(filesInTarR1)))
  readFilesTbl$fastqR1[readFilesTbl$tarFile==tarFile] <- filesInTarR1[idx]
  
  idx <- match( basename(readFilesTbl$PATH[readFilesTbl$tarFile==tarFile]),
                basename(dirname(filesInTarR2)))
  readFilesTbl$fastqR2[readFilesTbl$tarFile==tarFile] <- filesInTarR2[idx]
}


# Create a script to extract and rename files
outFilesR1 <- file.path(outDir,paste0(gsub("\\.","_",readFilesTbl$sampleID),"_R1.fastq.gz"))
outFilesR2 <- file.path(outDir,paste0(gsub("\\.","_",readFilesTbl$sampleID),"_R2.fastq.gz"))

script <- c(paste("tar -xf", readFilesTbl$tarFile, readFilesTbl$fastqR1, "-O >", outFilesR1),
            paste("tar -xf", readFilesTbl$tarFile, readFilesTbl$fastqR2, "-O >", outFilesR2))

source("processes/SLURMscript/createSLURMscript.R")

createSLURMscript(script,workdir = outDir,jobName = "extract")


# get RSEM result files

# StLa.T4.mix.1.genes.results
RSEMfiles <- file.path("/mnt/NOBACKUP/mariansc/share/RSEM",readFilesTbl$assembly,paste0(readFilesTbl$oldSampleID, ".genes.results"))
readFilesTbl$RSEMfilesDest <- paste0(gsub("\\.","_",readFilesTbl$sampleID), ".RSEMcounts.genes.txt")

file.copy(RSEMfiles,file.path(outDir,readFilesTbl$RSEMfilesDest))

# get trinity fasta files
spcs <- unique(readFilesTbl$assembly)
file.copy(file.path("/mnt/NOBACKUP/mariansc/share/trinity",spcs,paste0(spcs,".fasta")),
          file.path(outDir,paste0(spcs,".fasta")))
readFilesTbl$fastaDest <- paste0(readFilesTbl$assembly,".fasta")


# generate sample table
#
# Sample attributes:
  # Name: E.g. HoVu.W0.mix.1
  # Organism: E.g. "Hordeum vulgare"
  # Genotype: HoVu="variety Igri", BrDi="line BD1-1", StLa="wild type", MeNu1="wild type", NaSt="wild type"
  # "Single individual or mix": “mix”/”single”
  # (alternatively) Time after light on: “0 hours”/”8 hours”
  # Time: “0 (control)” / ”1 day” / ”4 weeks” / ”9 weeks”
  # Temperature: “17°C”/”6°C”
  # Day length: “12 hours”/”8 hours”
  # Light intensity: “150 µmol/m2s”/”50 µmol/m2s”

genotypes <- c(HoVu="variety Igri", BrDi="line BD1-1", StLa="wild type", MeNu1="wild type", NaSt="wild type")
organisms <- c(HoVu="Hordeum vulgare", BrDi="Brachypodium distachyon", StLa="Stipa lagascae", MeNu1="Melica nutans", NaSt="Nardus stricta")
timepoints <- c(`T-1`="0 (control)", T0="0 (control)", T1="1 day", T3="4 weeks",T4="9 weeks")

sampleTbl <- data.frame( sampleID = readFilesTbl$sampleID,
                         Organism = organisms[readFilesTbl$assembly],
                         Genotype = genotypes[readFilesTbl$assembly],
                         mix = ifelse(grepl("mix",readFilesTbl$sampleID),"mix","single"),
                         timeAfterLightOn = ifelse(grepl("W",readFilesTbl$sampleID),0, 8),
                         time = timepoints[readFilesTbl$TIMEPOINT],
                         temperature = ifelse(grepl("[DW]0",readFilesTbl$sampleID),17,6),
                         dayLength = ifelse(grepl("[DW]0",readFilesTbl$sampleID),12,8),
                         lightIntensity = ifelse(grepl("[DW]0",readFilesTbl$sampleID),"150 µmol/m2s","50 µmol/m2s"),
                         fastqR1 = paste0(gsub("\\.","_",readFilesTbl$sampleID),"_R1.fastq.gz"),
                         fastqR2 = paste0(gsub("\\.","_",readFilesTbl$sampleID),"_R2.fastq.gz"),
                         RSEMcounts = readFilesTbl$RSEMfilesDest,
                         assembly = readFilesTbl$fastaDest,
                         stringsAsFactors = F)

write.csv(sampleTbl,file = file.path(outDir,"sampleTbl.txt"),row.names = F)
                         
                                       
cat(sampleTbl$lightIntensity,sep="\n")

                         
                         

# setwd(outDir)
# system('md5sum *.gz > md5sum.txt')
# system('md5sum *.genes.txt > md5sum.genes.txt')
# system('md5sum *.fasta > md5sum.fasta.txt')
