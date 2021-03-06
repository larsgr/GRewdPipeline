source("R/fillTemplateFile.R")
source("processes/SLURMscript/createSLURMscript.R")

createTrimmoJob <- function( readFileFolders, outNames, outDir, jobArraySize, adapterFile){

  # stop if outDir already exists
  if(file.exists(outDir)){
    stop(paste0("Could not create job because directory ",outDir," already exists!"))
  }
  
  dir.create(outDir) # create output directory
  dir.create(file.path(outDir,"fastqc_output"))
  
  write.table( data.frame(outNames,readFileFolders),
               file = file.path(outDir,"sampleList.txt"),
               sep=" ", col.names=F, row.names=F, quote=F)
  
  file.copy("processes/trimmo/RUNtrimmo.sh", outDir, copy.mode=F)
  file.copy(adapterFile, outDir, copy.mode=F)
  
  script = fillTemplateFile( list(samplelist="sampleList.txt",arraysize=jobArraySize),
                            templateFile = "processes/trimmo/trimmo.Job.template.sh")
  createSLURMscript(script = script,workdir = normalizePath(outDir),jobName = "trimmo",
                    ntasks = 2, 
                    SBATCHoptions=paste0("#SBATCH -a 1-",jobArraySize))
}