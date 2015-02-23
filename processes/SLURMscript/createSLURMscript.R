source("R/fillTemplateFile.R")

# function that generates a SLURM script into which you can insert a script of your own

createSLURMscript <- function(script,workdir,jobName,ntasks=1,nodes=1,SBATCHoptions=""){
  txt <- fillTemplateFile(c(script=paste(script,collapse="\n"), workdir=workdir, ntasks=ntasks,
                            nodes=nodes,jobName=jobName,SBATCHoptions=SBATCHoptions),
                          templateFile="processes/SLURMscript/SLURMscript.template.sh")
  writeLines(txt,con = file.path(workdir,paste0(jobName,".job.sh")))  
}