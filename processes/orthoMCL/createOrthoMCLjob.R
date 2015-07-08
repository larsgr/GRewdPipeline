source("processes/SLURMscript/createSLURMscript.R")

createOrthoMCLjob <- function( outDir, 
                               proteomeFiles, 
                               taxon_codes, 
                               blastCPU=1,
                               blastArraySize=1){
  # stop if outDir already exists
  if(file.exists(outDir)){
    stop(paste0("Could not create job because directory ",outDir," already exists!"))
  }
  
  # check lengths
  if(length(proteomeFiles) != length(taxon_codes)){
    stop(paste0("Length of proteomeFiles and taxon_codes must be equal!"))
  }
  

  dir.create(outDir) # create output directory
  
  # normalize paths just in case
  outDir <- normalizePath(outDir)
  proteomeFiles <- normalizePath(proteomeFiles)
  
  # copy config file
  file.copy("processes/orthoMCL/orthomcl.config",outDir)
  # copy splitFasta script
  file.copy("processes/orthoMCL/splitFasta",outDir)
  
  prepScript <- paste( sep="\n",
                      "module load orthomcl",
                      "module load blast+",
                      "",
                      "mkdir compliantFasta",
                      "cd compliantFasta",
                      paste(sapply(seq_along(proteomeFiles),function(i){
                        paste("orthomclAdjustFasta", taxon_codes[i], proteomeFiles[i], "1")
                      }), collapse="\n"),
                      "",
                      "cd ..",
                      "cat compliantFasta/*.fasta > allProteomes.fasta",
                      "",
                      "makeblastdb -in allProteomes.fasta -dbtype prot -out allProteomes.db",
                      "",
                      "mkdir tmp",
                      paste("./splitFasta allProteomes.fasta",blastArraySize,"tmp/part") )
  
  
  createSLURMscript(script = prepScript, workdir = outDir, 
                    jobName = "prepOrthoMCL")
 
  
  blastScript <- paste(sep="\n",
                       paste("blastp",
                             "-query", sprintf("tmp/part_%04d.fa",1:blastArraySize),
                             "-db allProteomes.db",
                             "-outfmt 6",
                             "-out",sprintf("tmp/part_%04d.blastp",1:blastArraySize),
                             "-num_threads", blastCPU))
  write(blastScript, file.path(outDir,"commands.txt"))

  createSLURMarray(commandListFile = "commands.txt",
                   arraySize = blastArraySize,
                   workdir = outDir,
                   preScript="module load blast+",
                   jobName="blastOrthoMCL",
                   ntasks=blastCPU)
                           
  
  '
module load orthomcl
module load mcl

# concatenate the blast results
cat tmp/*.blastp > all_vs_all.out

# Database must be initialized first:
#orthomclInstallSchema orthomcl.config install_schema.log

echo "========= Step 8: orthomclBlastParser ========"

orthomclBlastParser all_vs_all.out compliantFasta >> similarSequences.txt

echo "DONE Blast output parsed to orthoMCL format"
  
  
echo "========= Step 9: orthomclLoadBlast  ==========="
  
orthomclLoadBlast orthomcl.config  similarSequences.txt

echo "DONE similar sequences loaded into DB"


echo "========= Step 10: orthomclPairs ========="

orthomclPairs orthomcl.config orthomclPairs.log cleanup=yes

echo "DONE orthomclPairs"


echo "========== Step 11: orthomclDumpPairsFiles ========"

orthomclDumpPairsFiles orthomcl.config

echo "DONE orthomclDumpPairsFiles"


echo "========== Step 12: mcl ========"

mcl mclInput --abc -I 1.5 -o mclOutput

echo "DONE mcl"


echo "========== Step 13: orthomclMclToGroups =========="

orthomclMclToGroups grp 100000 < mclOutput > groups.txt

echo "DONE orthomclMclToGroups"
  ' -> orthoMCLscript
  
  createSLURMscript(script = orthoMCLscript, workdir = outDir,
                    jobName = "orthoMCL")
  
}
