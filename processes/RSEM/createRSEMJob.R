source("processes/SLURMscript/createSLURMscript.R")

createRSEMJob <- function( outDir, leftReadFiles, rightReadFiles, outputPrefixes, 
                           transcriptsFile, jobName = "RSEM", CPU=10, arraySize=1,
                           prepCPU=1, seqType="fq", SS_lib_type="RF"){
  # stop if outDir already exists
  if(file.exists(outDir)){
    stop(paste0("Could not create job because directory ",outDir," already exists!"))
  }
  
  nSamples <- length(outputPrefixes)
  if( nSamples < 1 | length(leftReadFiles) != nSamples | length(rightReadFiles) != nSamples ){
    stop("Length of leftReadFiles, rightReadFiles and outputPrefixes must be equal and > 0")
  }
  
  if(nSamples < arraySize){
    stop("arraySize can not be greater than number of samples")
  }
  
  dir.create(outDir) # create output directory
  
  # normalize paths just in case
  outDir <- normalizePath(outDir)
  transcriptsFile <- normalizePath(transcriptsFile)
  
  
  
  # load modules:  
  preScript <- paste( sep="\n",
                      "module load trinity",
                      "module load rsem",
                      "module load samtools",
                      "")

  # create transcriptome prep script
  prepScript <- paste( sep="\n", preScript,
                       paste("/local/genome/packages/trinity/2.0.2/util/align_and_estimate_abundance.pl",
                             "--transcripts", transcriptsFile,
                             "--est_method RSEM",
                             "--aln_method bowtie",
                             "--trinity_mode",
                             "--prep_reference",
                             "--thread_count", prepCPU))
  
  createSLURMscript(script = prepScript, workdir = outDir, 
                    jobName = paste0("prep",jobName),ntasks = prepCPU)
  
  # Build command for each sample:
  commandList <- character(0)
  for( i in 1:nSamples){
    cmd <- paste(
      "/local/genome/packages/trinity/2.0.2/util/align_and_estimate_abundance.pl ",
      "--transcripts", transcriptsFile,
      "--seqType", seqType,
      "--left", leftReadFiles[i],
      "--right", rightReadFiles[i],
      "--est_method RSEM", 
      "--aln_method bowtie",
      "--SS_lib_type", SS_lib_type,
      "--thread_count", CPU,
      "--output_dir", outDir,
      "--trinity_mode",
      "--output_prefix", outputPrefixes[i])
    
    # add command to script
    commandList <- c(commandList,cmd)
  }
  
  
  # write list of commands to file
  writeLines(commandList,con=file.path(outDir,"commandList.txt"))
  
  # wrap into a slurm job script and write to file
  createSLURMarray(commandListFile="commandList.txt",
                   arraySize=arraySize,
                   workdir=outDir,
                   preScript=preScript,
                   jobName=jobName,
                   ntasks=CPU)
}
