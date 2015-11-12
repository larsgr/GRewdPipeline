#
# NOTE: system("source /etc/profile.d/modules.sh && sbatch jobscript.sh") does not work!
# 
# The module command is not available to "jobscript.sh".. why not?
# apparently the job script does inherit the environment where sbatch is executed.
#
# Need to add #!/bin/bash –l at the start of each job script...

submitJob <- function(jobScript, deps=NULL){
  if(is.null(deps)){
    cmd <- paste("sbatch",jobScript)
  } else {
    cmd <- paste0("sbatch -d afterok:",paste(deps,collapse=":")," ",jobScript)    
  }
  cat("Submit job command:",cmd,"\n")
  jobIDstr <- system(paste("source /etc/profile.d/modules.sh && module load slurm &&",cmd),intern=T)
  cat(jobIDstr,"\n")
  # Output from sbatch: "Submitted batch job 1284123"
  if( !( length(jobIDstr) == 1 & grepl("^Submitted batch job [0-9]+$",jobIDstr))){
    stop("Job submit failed!")
  } 
  jobID <- stringr::str_extract(jobIDstr,"[0-9]+")  
  return(jobID)
}

# 
# orthoMCL.prepOrthoMCL.jobID <- submitJob(
#   jobScript = "/mnt/NOBACKUP/mariansc/share/orthos/orthoMCL/prepOrthoMCL.job.sh",
#   deps=NULL)
# 
# orthoMCL.blastOrthoMCL.jobID <- submitJob(
#   jobScript = "/mnt/NOBACKUP/mariansc/share/orthos/orthoMCL/blastOrthoMCL.job.sh",
#   deps=orthoMCL.prepOrthoMCL.jobID)
# 
# 
# orthoMCL.orthoMCL.jobID <- submitJob(
#   jobScript = "/mnt/NOBACKUP/mariansc/share/orthos/orthoMCL/orthoMCL.job.sh",
#   deps=orthoMCL.blastOrthoMCL.jobID)
# 
# 
# ExTbls.jobID <- submitJob(
#   jobScript = "/mnt/NOBACKUP/mariansc/share/orthos/exprTbls/ExTbls.job.sh",
#   deps=orthoMCL.orthoMCL.jobID)
# 
# runDESeq.jobID <- submitJob(
#   jobScript = "/mnt/NOBACKUP/mariansc/share/orthos/DESeq/runDESeq.job.sh",
#   deps=ExTbls.jobID)
# 
# grpFastas.RJob.jobID <- submitJob(
#   jobScript = "/mnt/NOBACKUP/mariansc/share/orthos/grpFastas/RJob.job.sh",
#   deps=orthoMCL.orthoMCL.jobID)
# 
# grpCDSFastas.jobID <- submitJob(
#   jobScript = "/mnt/NOBACKUP/mariansc/share/orthos/grpCDSFastas/grpCDSFastas.job.sh",
#   deps=orthoMCL.orthoMCL.jobID)
# 
# prepMAFFT.jobID <- submitJob(
#   jobScript = "/mnt/NOBACKUP/mariansc/share/orthos/grpAligned/prepMAFFT.job.sh",
#   deps=grpFastas.RJob.jobID)
# 
# MAFFT.jobID <- submitJob(
#   jobScript = "/mnt/NOBACKUP/mariansc/share/orthos/grpAligned/MAFFT.job.sh",
#   deps=prepMAFFT.jobID)
# 
# pal2nal.jobID <- submitJob(
#   jobScript = "/mnt/NOBACKUP/mariansc/share/orthos/pal2nal/pal2nal.job.sh",
#   deps=c(MAFFT.jobID,grpCDSFastas.jobID))
# 
# treesNuc.jobID <- submitJob(
#   jobScript = "/mnt/NOBACKUP/mariansc/share/orthos/treesNuc/nucTree.job.sh",
#   deps=pal2nal.jobID)

splitGroups.jobID <- submitJob(
  jobScript = "/mnt/NOBACKUP/mariansc/share/orthos/splitGroups/splitGroups.job.sh",
  deps=treesNuc.jobID)

PAML.jobID <- submitJob(
  jobScript = "/mnt/NOBACKUP/mariansc/share/orthos/PAML/codeml.job.sh",
  deps=splitGroups.jobID)
