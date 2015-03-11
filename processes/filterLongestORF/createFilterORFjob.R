source("processes/SLURMscript/createSLURMscript.R")

createFilterORFjob <- function( outDir, ORFfiles, outPrefix, jobName = "filterORFs",
                                arraySize = length(ORFfiles) ){
  
  # stop if outDir already exists
  if(file.exists(outDir)){
    stop(paste0("Could not create job because directory ",outDir," already exists!"))
  }
  
  nFiles <- length(ORFfiles)
  if(nFiles < arraySize){
    stop("arraySize can not be greater than number of files")
  }
  
  dir.create(outDir) # create output directory
  
  file.copy(from = "processes/filterLongestORF/filterLongestORF.R",outDir)
  
  # Build command for each file:
  commandList <- character(0)
  for( i in 1:nFiles){
    cmd <- paste(
      "Rscript filterLongestORF.R",
      ORFfiles[i],
      paste0(outPrefix[i],".longest.fasta"),
      paste0(outPrefix[i],".tbl") )
    
    # add command to script
    commandList <- c(commandList,cmd)
  }
  
  # write list of commands to file
  writeLines(commandList,con=file.path(outDir,"commandList.txt"))
   
  
  # wrap into a slurm job script and write to file
  createSLURMarray(commandListFile="commandList.txt",
                   arraySize=arraySize,
                   workdir=normalizePath(outDir),
                   preScript="module load R\n",
                   jobName=jobName )
  
  
}
