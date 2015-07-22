source("processes/SLURMscript/createSLURMscript.R")
source("processes/RJob/RJob.R")


MAFFTJob <- function( inFastaDir,
                      outDir,
                      jobName = "MAFFT",
                      arraySize = 100,
                      nThreads = 1 ){
  list(
    outDir = outDir,
    jobName = jobName,
    params = list(
      nThreads = nThreads,
      arraySize = arraySize) ,
    input = list(
      inFastaDir = inFastaDir )
  ) -> job
  class(job) <- append(class(job),"MAFFTJob")
  return(job)
}

generateScript.MAFFTJob = function( job, overwrite=FALSE ){
  with(job, {
    
    # make a script "prepMAFFT.job.sh" that generates the commands.txt
    generateScript(
      RJob(data=list(inFastaDir=input$inFastaDir),outDir=outDir,jobName="prepMAFFT",
           FUN=function(){
             fastaFiles <- dir(data$inFastaDir,pattern="\\.fasta$",full.names = T)
             alnFiles <- sub("fasta$","aln",basename(fastaFiles))
             idx <- !file.exists(alnFiles) # skip already aligned files
             writeLines(paste("mafft --quiet",fastaFiles[idx],">",alnFiles[idx]),"commands.txt")
           } ),
      overwrite = overwrite
      )
    
    
    # generate SLURM job script
    createSLURMarray(commandListFile = "commands.txt",
                     preScript = "module load mafft/7.130",
                     workdir = outDir,
                     jobName = jobName,
                     arraySize = params$arraySize,
                     ntasks = params$nThreads)
    
  }) 
}

generateScript <- function(job, ...) UseMethod( "generateScript" )
# generateScript <- function(job){
#   UseMethod( "generateScript", job )
# }
