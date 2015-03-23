source("R/fillTemplateFile.R")

# function for converting list of parameters to SBATCH parameters
makeSbatchOptions <- function (params, SBATCHoptions="") {
  for( i in seq_along(params)){
    paramName <- gsub("_","-",names(params)[i])
    paramValue <- params[[i]]
    if(is.logical(paramValue)){
      if(paramValue==TRUE)
        SBATCHoptions <- paste(
          paste0("#SBATCH --",paramName),
          SBATCHoptions, sep="\n")
    } else {
      SBATCHoptions <- paste(
        paste0("#SBATCH --",paramName,"=",paramValue),
        SBATCHoptions, sep="\n")
    }
  }
  return(SBATCHoptions)
}

# function that generates a SLURM script into which you can insert a script of your own


createSLURMscript <- function(script,workdir,jobName,ntasks=1,nodes=1,SBATCHoptions="",...){
  
  # pass additional ... parameters to sbatch
  SBATCHoptions <- makeSbatchOptions(params = list(...),SBATCHoptions)
  
  
  txt <- fillTemplateFile(c(script=paste(script,collapse="\n"),
                            workdir=normalizePath(workdir),
                            ntasks=ntasks,
                            nodes=nodes,
                            jobName=jobName,
                            SBATCHoptions=SBATCHoptions),
                          templateFile="processes/SLURMscript/SLURMscript.template.sh")
  writeLines(txt,con = file.path(workdir,paste0(jobName,".job.sh")))  
}



# create a SLURM array script that executes a list of commands
createSLURMarray <- function(commandListFile,arraySize,workdir,preScript="",
                             jobName="array",ntasks=1,nodes=1,...){
  
  
  # pass additional ... parameters to sbatch
  SBATCHoptions <- makeSbatchOptions(params = list(...))
  
  
  txt <- fillTemplateFile(c(preScript=paste(preScript,collapse="\n"), 
                            workdir=normalizePath(workdir), 
                            ntasks=ntasks, 
                            arraySize=arraySize,
                            commandList=commandListFile,
                            nodes=nodes,
                            jobName=jobName,
                            SBATCHoptions=SBATCHoptions),
                          templateFile="processes/SLURMscript/SLURMarray.template.sh")
  writeLines(txt,con = file.path(workdir,paste0(jobName,".job.sh")))  
}
