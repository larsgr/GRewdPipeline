loadOrthoGrpsTable <- function(orthoGrpFile){
  # Load file as a character vector.
  ortVec <- scan(file=orthoGrpFile,what=character(0))  

  # This vector is a mix of group IDs and gene IDs, which has to be seperated.  
  
  # Detect which are gene references
  isGeneId <- grepl("\\|",ortVec) # (genes have '|' in the name, group IDs don't)
  whichIsGrpId <- which(!isGeneId) # if its not a gene then its a group
  grpIds <- ortVec[whichIsGrpId] # store group IDs in seperate vector
  
  # get number of genes per group
  genesPerGrp <- (c(whichIsGrpId[-1],length(ortVec)+1) - whichIsGrpId) - 1

  grpVec <- rep(sub(":","",grpIds),times=genesPerGrp)
  
  # store gene IDs in seperate vector
  genVec <- ortVec[isGeneId]
  
  grpTable <- data.frame(seqID=genVec,grpID=grpVec,stringsAsFactors = F)
  
  return(grpTable)
}

OrthoGrpTableToArray <- function(grpTable){
  genSpcArray <- stringr::str_split_fixed(grpTable$seqID,"\\|",2)
  
  grps <- tapply(genSpcArray[ ,2],list(grp=grpTable$grpID,spc=genSpcArray[ ,1]), function(x){x})
  return(grps)
}

loadOrthoGrpsArray <- function(orthoGrpFile){
  OrthoGrpTableToArray( loadOrthoGrpsTable( orthoGrpFile ) )
}
  

grpToDF <- function(grp){
  maxParalogs <- max(sapply(grp,length))
  as.data.frame(lapply(grp,function(seq){
    c(seq,rep(NA_character_,maxParalogs-length(seq)))
  }),stringsAsFactors=F)
}

grpToChar <- function(grp){
  seqsPerSpecies <- sapply(grp,length)
  paste(rep(names(seqsPerSpecies),seqsPerSpecies),unlist(grp),sep="|")
}