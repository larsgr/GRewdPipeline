source("processes/SLURMscript/createSLURMscript.R")

createExprTablesJob <- function(outDir,orthoGrpFile,RSEMout,jobName="ExTbls"){
  # stop if outDir already exists
  if(file.exists(outDir)){
    stop(paste0("Could not create job because directory ",outDir," already exists!"))
  }
  
  dir.create(outDir)
  
  file.copy("processes/makeExprTables/makeExprTables.R",outDir)
  
  save(outDir,orthoGrpFile,RSEMout,file = file.path(outDir,"params.RData"))
  
  script <- "
module load R

Rscript makeExprTables.R
"
  
  createSLURMscript(script,workdir=outDir,jobName=jobName)
}