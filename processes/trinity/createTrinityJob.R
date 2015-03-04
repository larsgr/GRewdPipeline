source("R/fillTemplateFile.R")
source("processes/SLURMscript/createSLURMscript.R")

createTrinityJob <- function( leftReadFiles, rightReadFiles, outDir, trinityOutputName,
                              jobName = "trinity", seqType="fq", SS_lib_type="RF", 
                              CPU=10, max_memory="calculate"){
  
  # stop if outDir already exists
  if(file.exists(outDir)){
    stop(paste0("Could not create job because directory ",outDir," already exists!"))
  }
  
  dir.create(outDir) # create output directory
  
  if(max_memory=="calculate"){
    # require 4x fq file size
    memUsageGB <- max(2,4*sum(file.info(leftReadFiles)$size) %/% 1000000000)
    max_memory <- paste0(memUsageGB,"G")
  }
  

  script <- paste( sep="\n",
                   
                   "module load trinity",
                   "module load bowtie",
                   "module load samtools",
                   "",
                   paste(
                     "Trinity",
                     "--seqType", seqType,
                     "--SS_lib_type", SS_lib_type,
                     "--left", paste(leftReadFiles,collapse=','),
                     "--right", paste(rightReadFiles,collapse=','),
                     "--CPU", CPU,
                     "--max_memory", max_memory,
                     "--full_cleanup" ),
                   "",
                   paste0('
if [ -f trinity_out_dir.Trinity.fasta ];
then
  echo "CMD: mv trinity_out_dir.Trinity.fasta ',trinityOutputName,'"
  mv trinity_out_dir.Trinity.fasta ',trinityOutputName,'
else
  echo "Trinity assembly file does not exists. Check log for errors and try again"
  echo "[$(date +"%Y-%m-%d %H:%M:%S")] Job finished with errors" >> trinity.$SLURM_JOB_ID.started
  mv ',jobName,'.$SLURM_JOB_ID.started ',jobName,'.$SLURM_JOB_ID.failed
  exit 1
fi')
                   )
  
  
  createSLURMscript(script = script,workdir = normalizePath(outDir),jobName = jobName,
                    ntasks = CPU, partition="hugemem", mem=max_memory)
}
