# Convert RSEM output files to more handy expression tables
#
# For each species:
#   Create tables containing only genes that are included in an ortholog group
#     Change the names of the genes so that they match those in the groups
#   Store the tables
#
# Also create cross species tables
#   singleton genes (1:1 orthologs) table
#   sum of paralogs table
#
# All tables should have one version with raw read counts and one with FPKM




source("/mnt/users/lagr/networkSimilarity/R/loadOrthoGroups.R")


load("params.RData")

# params.RData contains:
#   outDir: output directory
#   orthoGrpFile: e.g."/mnt/NOBACKUP/mariansc/share/orthos/orthoMCL/groups.txt"
#   RSEMout as defined in pipeline.R

grps <- loadOrthoGrpsArray(orthoGrpFile = orthoGrpFile)


lapply( RSEMout, function(files.RSEM.spc){
  lapply( files.RSEM.spc, function(files.RSEM.spc.T){
    cat("Reading",files.RSEM.spc.T,"\n")
    read.table(files.RSEM.spc.T, header=T)
  })  
}) -> RSEMexpr


spcs <- names(RSEMexpr)
spcs <- setNames(spcs,spcs)

# get number of genes per species per group
lens <- unlist(lapply(grps,length))
dim(lens) <- dim(grps)
dimnames(lens) <- dimnames(grps)


# get singleton groups
grps11 <- apply(lens==1,1,all)


# create cross species tables
# singleton genes (1:1 orthologs) table
grpSingletonTbls <- list()
# Sum of paralogs per group
grpSumTbls <- list()

# make one table for FPKM and one with raw counts
countTypes <- c("FPKM","expected_count")

myWriteTbl <- function(x,fileName){
  cat("Writing",fileName,"\n")  
  write.table(x,file = fileName,quote=F,sep="\t")
}

for( spc in spcs ){
  gene_ids <- gsub("|","_",RSEMexpr[[spc]][[1]]$gene_id, fixed = T)
  gene_ids_grps <- unlist(grps[,spc])
  isInGrps <- gene_ids %in% gene_ids_grps
  
  for( countType in countTypes){
    countTbl <- sapply(RSEMexpr[[spc]],function(X){X[[countType]][isInGrps]})
    rownames(countTbl) <- gene_ids[isInGrps]
    
    myWriteTbl(countTbl, file.path(outDir,paste0(spc,"_",countType,"Tbl.txt")) )
    
    # get singleton genes:
    gene_idx_singleton <- match(unlist(grps[grps11,spc]),rownames(countTbl))
    singletonCountTbl <- countTbl[gene_idx_singleton, ]
    rownames(singletonCountTbl) <- rownames(grps)[grps11]
    grpSingletonTbls[[countType]] <- cbind(grpSingletonTbls[[countType]],singletonCountTbl)
    
    # calculate paralog expression sums
    sumCountTbl <- matrix(0,nrow = nrow(grps), ncol = ncol(countTbl),
                          dimnames = list(rownames(grps),colnames(countTbl)))
    hasOrtho <- lens[,spc]>0
    sumCountTbl[hasOrtho,] <- t(sapply(grps[hasOrtho,spc],function(id){
      colSums(countTbl[id, ,drop=F]) 
    }))
    
    grpSumTbls[[countType]] <- cbind(grpSumTbls[[countType]],sumCountTbl)
    
  }  
  
}

myWriteTbl(grpSingletonTbls$FPKM, file.path(outDir,"grpSingletonFPKMTbl.txt"))
myWriteTbl(grpSingletonTbls$expected_count, file.path(outDir,"grpSingletonCountTbl.txt"))
myWriteTbl(grpSumTbls$FPKM, file.path(outDir,"grpSumFPKMTbl.txt"))
myWriteTbl(grpSumTbls$expected_count,file.path(outDir,"grpSumCountTbl.txt"))
