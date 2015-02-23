source("pipeline/R/fillTemplateFile.R")

# function that generates a SLURM script into which you can insert a script of your own

createSLURMscript <- function(script,workdir,jobName,ntasks=1,nodes=1){
  txt <- fillTemplateFile(c(script=script, workdir=workdir, ntasks=ntasks,
                            nodes=nodes, output=output,jobName=jobName),
                          templateFile="pipeline/processes/SLURMscript/SLURMscript.template.sh")
  writeLines(txt,con = file.path(workdir,paste0(jobName,".job.sh")))  
}