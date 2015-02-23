source("R/fillTemplateFile.R")
source("processes/SLURMscript/createSLURMscript.R")

createTrimmoJob <- function( readFileFolders, outNames, outDir, jobArraySize, adapterFile){
  
  dir.create(outDir) # create output directory
  dir.create(file.path(outDir,"fastqc_output"))
  
  write.table( data.frame(outNames,readFileFolders, rep(".",length(readFileFolders))),
               file = file.path(outDir,"sampleList.txt"),
               sep=" ", col.names=F, row.names=F, quote=F)
  
  file.copy("processes/trimmo/RUNtrimmo.sh", outDir)
  file.copy(adapterFile, outDir)
  
  script = fillTemplateFile( list(samplelist="sampleList.txt",arraysize=jobArraySize),
                            templateFile = "processes/trimmo/trimmo.Job.template.sh")
  createSLURMscript(script = script,workdir = normalizePath(outDir),jobName = "trimmo",
                    ntasks = 2, 
                    SBATCHoptions=paste0("#SBATCH -a 1-",jobArraySize))
}