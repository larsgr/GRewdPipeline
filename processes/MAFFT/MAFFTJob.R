source("processes/SLURMscript/createSLURMscript.R")


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

generateScript.MAFFTJob = function( job ){
  with(job, {
    # stop if outDir already exists
    if(file.exists(outDir)){
      stop(paste0("Could not create job because directory ",outDir," already exists!"))
    }
    
    dir.create(outDir, recursive = T) # create output directory
    
    #### WARNING: The preScript is executed in paralell, but the commands.txt should only
    ####          be generated once!
    ####
    #### TODO: Should set quiet parameter and threads parameter
    # generate the commands in the preScript
    # The command is:
    #   mafft $1 > $2.aln
    # where $1 is the full filename of the input fasta and $2 is the basename 
    # without .fasta extension
    preScript = paste(sep="\n",
                      "module load mafft/7.130",
                      "",
                      paste("find", input$inFastaDir, "-iname '*.fasta'",
                            "|",
                            "perl -pe's/(.*\\/([^\\/]*)(\\.fasta))/mafft $1 > $2.aln/'",
                            ">",
                            "commands.txt") )
                      
    
    
    # generate SLURM job script
    createSLURMarray(commandListFile = "commands.txt",
                     preScript = preScript,
                     workdir = outDir,
                     jobName = jobName,
                     arraySize = params$arraySize,
                     ntasks = params$nThreads)
    
  }) 
}

generateScript <- function(job){
  UseMethod( "generateScript", job )
}
