source("pipeline/R/fillTemplateFile.R")
source("pipeline/processes/SLURMscript/createSLURMscript.R")

# Create an R script that sources an R file then calls a function using the 
# commandline arguments as arguments for the function call.
createRscriptWrapper <- function(sourceFile, funName, outFile){
  txt <- fillTemplateFile(c(RSourceFile=sourceFile, FunctionName=funName),
                   templateFile = "pipeline/processes/RScript/wrapper.template.R")
  writeLines(txt,con = outFile)
}


# TODO: make a function that that generates Rscript launcher that passes the arguments
createRscriptJob <- function(outDir,sourceFile, funName, args, jobName=funName){

  wrapperName <- paste0(funName,".wrapper.R")
  createRscriptWrapper(sourceFile, funName, 
                       outFile=file.path(outDir,wrapperName))
  script <- sprintf("
module load R

Rscript %s %s
",wrapperName, paste(args,collapse=" "))
  
  createSLURMscript(script,workdir=outDir,jobName=jobName)
}


# example?:
#
# createRscriptJob(outDir = normalizePath("~/GRewd/pipelineData/orthoGrpFastas"),
#                  sourceFile = "~/GRewd/pipeline/processes/getOrthoGrpFastas/getOrthoGrpFastas.R",
#                  funName = "getOrthoGrpFastas", 
#                  args=list(grpFile="~/GRewd/orthoMCL/groups.txt",
#                            fastaFile="/mnt/users/mariansc/my_orthomcl_dir/allProteomes.fasta",
#                            outDir="~/GRewd/pipelineData/orthoGrpFastas/fastas"))
