source("processes/SLURMscript/createSLURMscript.R")

RJob <- function( data=NULL, FUN, outDir, jobName = "RJob", nThreads=1){
  if( !is.function(FUN) ) stop("FUN must be a function!")
  
  list(
    outDir = outDir,
    jobName = jobName,
    params = list(
      nThreads = nThreads,
      data = data,
      FUN = FUN )
  ) -> job
  class(job) <- c("RJob","list")
  return(job)
}

generateScript.RJob = function( job, overwrite=FALSE ){
  with(job, {
    # stop if outDir already exists
    if(file.exists(outDir) & (overwrite==FALSE) ){
      stop(paste0("Could not create job because directory ",outDir," already exists!"))
    }
    
    dir.create(outDir, recursive = T) # create output directory
    
    scriptPath <- file.path(outDir,"script")
    dir.create( scriptPath ) # create script directory
    
    # copy runFUN.R script
    file.copy("processes/RJob/runRFUN.R", scriptPath)
    
    # save the data
    saveRDS(params$data, file.path(scriptPath,"data.RDS"))
    saveRDS(params$FUN, file.path(scriptPath,"FUN.RDS"))
    
    # generate the script:
    script <- paste( sep="\n",
                     "module load R",
                     "",
                     "Rscript script/runRFUN.R" )
    
    # generate SLURM job script
    createSLURMscript(script = script,
                     workdir = outDir,
                     jobName = jobName,
                     ntasks = params$nThreads)
  }) 
}

generateScript <- function(job, ...) UseMethod( "generateScript" )

