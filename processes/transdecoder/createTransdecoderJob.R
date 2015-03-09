source("processes/SLURMscript/createSLURMscript.R")

createTransdecoderJob <- function( outDir, transcriptsFile,
                              jobName = "transdecoder", strandSpecific=TRUE,
                              minProtLength = 30,
                              CPU=10){
  
  # stop if outDir already exists
  if(file.exists(outDir)){
    stop(paste0("Could not create job because directory ",outDir," already exists!"))
  }
  
  dir.create(outDir) # create output directory
  
  hmmFile <- "/local/genome/packages/transdecoder/rel16JAN2014/pfam/Pfam-AB.hmm.bin"
  
  script <- paste( sep="\n",
                   
                   "module load transdecoder/rel16JAN2014",
                   "",
                   paste(
                     "TransDecoder",
                     "-t", normalizePath(transcriptsFile),
                     "--search_pfam", hmmFile,
                     "-m", minProtLength,
                     ifelse(strandSpecific,"-S",""),
                     "--CPU", CPU
                     )
                  )
  
  
  createSLURMscript(script = script,workdir = normalizePath(outDir),jobName = jobName,
                    ntasks = CPU)
}
