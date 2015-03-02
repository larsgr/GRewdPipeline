source("R/fillTemplateFile.R")
source("processes/SLURMscript/createSLURMscript.R")

createTrinityJob <- function( leftReadFiles, rightReadFiles, outDir,
                              seqType="fq", SS_lib_type="RF", CPU=10, max_memory="calculate"){
  
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
                     "--full_cleanup" ))
  
  
  createSLURMscript(script = script,workdir = normalizePath(outDir),jobName = "trinity",
                    ntasks = CPU, partition="hugemem", mem=max_memory)
}
