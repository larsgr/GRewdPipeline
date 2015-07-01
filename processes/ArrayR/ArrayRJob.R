source("processes/SLURMscript/createSLURMscript.R")

ArrayRJob <- function( x, FUN, outDir, jobName = "ArrayRJob", commonData=NULL ,
                       nThreads=1, arraySize=100 ){
  if( !is.function(FUN) ) stop("FUN must be a function!")
  if( length(x) < 1 ) stop("x must have length > 0!")

  list(
    outDir = outDir,
    jobName = jobName,
    params = list(
      nThreads = nThreads,
      x = as.list(x),
      commonData = commonData,
      FUN = FUN,
      arraySize = arraySize)
  ) -> job
  class(job) <- c("ArrayRJob","list")
  return(job)
}

generateScript.ArrayRJob = function( job ){
  with(job, {
    # stop if outDir already exists
    if(file.exists(outDir)){
      stop(paste0("Could not create job because directory ",outDir," already exists!"))
    }
    
    dir.create(outDir, recursive = T) # create output directory

    scriptPath <- file.path(outDir,"script")
    dir.create( scriptPath ) # create script directory
    
    # copy runFUN.R script
    file.copy("processes/ArrayR/runFUN.R", scriptPath)
    
    # save the data
    saveRDS(params$x, file.path(scriptPath,"x.RDS"))
    saveRDS(params$commonData, file.path(scriptPath,"commonData.RDS"))
    saveRDS(params$FUN, file.path(scriptPath,"FUN.RDS"))
    
    # generate the commands:
    commandsList <- paste("Rscript script/runFUN.R", 1:length(params$x))
    
    # write the commands to file
    writeLines(commandsList, file.path(outDir,"script/commands.txt"))
    
    # generate SLURM job script
    createSLURMarray(commandListFile = "script/commands.txt",
                     preScript = "module load R",
                     workdir = outDir,
                     jobName = jobName,
                     arraySize = params$arraySize,
                     ntasks = params$nThreads)
  }) 
}

generateScript <- function(job){
  UseMethod( "generateScript", job )
}

# 
# arrayRJob <- ArrayRJob( 1:3, function(x){
#   cat(" x =",x,"\n")
#   write(paste(commonData$msg,x),file = paste0(LETTERS[x],".txt"))
# }, outDir="~/temp/testArrayRJob", commonData = list(msg="Hello world!"))
# 
# generateScript(arrayRJob)
