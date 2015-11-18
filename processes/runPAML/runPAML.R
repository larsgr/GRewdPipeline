library(ape)
library(phangorn) # getClans()

RLinuxModules::moduleInit()
RLinuxModules::module("load paml/4.7a")


#############
#
# Helper functions:
#
isClan <- function(tree, isInGrp){
  clans = getClans(tree)
  any(apply(clans,1,function(clan){
    all(clan==isInGrp)
  }))
}

readFasta2phylip <- function (alnFile, seqIDs.splitClan) {
  # read fasta file
  txt.fasta <- readLines(alnFile)
  # convert to into named character vector
  headerLines <- grepl("^>",txt.fasta)
  seqIDs <- sub("^>","",txt.fasta[headerLines]) # OBS assumes only seqID on header
  seqs <- tapply(txt.fasta[!headerLines],
                 seqIDs[cumsum(headerLines)[!headerLines]],
                 paste,collapse="")
  
  # count bases
  nBases <- nchar(seqs[1])
  # if( !all( sapply(seqs,nchar) == nBases) ) # check if all sequences have same length
  # if( (nBases %% 3) != 0 ) # check if nBases is factor of 3
  
  # if grp has been split, remove the uneeded sequences
  seqs <- seqs[seqIDs.splitClan]
  
  txt.phylip <- c(paste(length(seqs),nBases),
                  paste(names(seqs),seqs,sep="\n"))
  
  return(txt.phylip)
}

ctlFileTxt <- function(seqfile,
                       treefile,  
                       outfile,
                       
                       noisy = 0,
                       verbose = 0,
                       runmode = 0,
                       
                       seqtype = 1,  # 1:codons; 2:AAs; 3:codons-->AAs
                       CodonFreq = 2,  # 0:1/61 each, 1:F1X4, 2:F3X4, 3:codon table
                       clock = 0,   # 0:no clock, 1:global clock; 2:local clock; 3:TipDate
                       aaDist = 0,  # 0:equal, +:geometric; -:linear, 1-6:G1974,Miyata,c,p,v,a
                       model = 2,
                       # models for codons:
                       # 0:one, 1:b, 2:2 or more dN/dS ratios for branches
                       
                       NSsites = 2,  # 0:one w;1:neutral;2:selection; 3:discrete;4:freqs;
                       # 5:gamma;6:2gamma;7:beta;8:beta&w;9:beta&gamma;
                       # 10:beta&gamma+1; 11:beta&normal>1; 12:0&2normal>1;
                       # 13:3normal>0
                       
                       icode = 0,  # 0:universal code; 1:mammalian mt; 2-10:see below
                       fix_kappa = 0,  # 1: kappa fixed, 0: kappa to be estimated
                       kappa = 3,  # initial or fixed kappa
                       fix_omega,  # 1: omega or omega_1 fixed, 0: estimate 
                       omega = 1,  # initial or fixed omega, for codons or codon-based AAs
                       
                       fix_alpha = 1,  # 0: estimate gamma shape parameter; 1: fix it at alpha
                       alpha = 0,  # initial or fixed alpha, 0:infinity (constant rate)
                       Malpha = 0,  # different alphas for genes
                       ncatG = 8,  # # of categories in dG of NSsites models
                       
                       getSE = 0,  # 0: don't want them, 1: want S.E.s of estimates
                       RateAncestor = 0,  # (0,1,2): rates (alpha>0) or ancestral states (1 or 2)
                       
                       Small_Diff = .5e-6,
                       cleandata = 0,
                       method = 0   # 0: simultaneous; 1: one branch at a time
                       )
{
  params <- as.list(environment())
   
  return( paste( names(params), "=", params,  collapse="\n") )
}

# mark the most recent common ancestor of tip matching pattern
markTree <- function(tree, pattern){
  tipIdx <- grep(pattern,tree$tip.label)
  if( length(tipIdx) == 1){
    # leaf node (add mark to tip.label)
    tree$tip.label[tipIdx] <- paste0(tree$tip.label[tipIdx],"#1")
  } else {
    # internal node
    idx <- getMRCA(tree,tipIdx)
    tree$node.label[idx-length(tree$tip.label)] <- "#1"    
  }
  return(tree)
}

# warning assumes the existance of the file "aln.phylip"
runCodeML <- function( markedTree, outFile, ... ){
  if(file.exists(outFile)){
    cat("Skipping codeml for existing file",outFile,"\n")    
  } else{
    outTmp <-  paste0(outFile,".tmp")
    write.tree(markedTree, file = "tree")
    write( ctlFileTxt(seqfile = "aln.phylip",
                      treefile = "tree",
                      outfile = outTmp,
                      ...),
           file = "codeml.ctl")  
    
    cat("Running codeml for",outFile,"\n")

    if(file.exists(outTmp)){
      unlink(outTmp) # remove old tmp file
    }
    
    system("codeml")

    if(file.exists(outTmp)){
      file.rename(outTmp,outFile)
    }
  }
  
}



# For each valid tree:
#   Translate alignment to phylip format
#     add a line with two numbers (number of seqs and number of aligned nucleotides) on the start of the file
#     Remove “>” infront of names
#   Remove BS numbers
#   for each hypothesesis: 
#     mark the branches that shall be tested and save tree file
#     Generate ctl file pointing to correct tree file and with correct hypothesis
#     Run codeml



# for each a valid tree:
# Run codeml based on the different hypothesis
# Store results under a folder for each hypothesis
#   name of result file is grp1234.out
# put temporary files in working directory
doPAML <- function(tree, grpID, codonAlnPath, outGrpDir){
  # Find and load the corresponding alignment fasta file
  # Split group IDs have a ".XX" suffix, get the base group ID:
  grpID.noSplit <- sub("\\.[0-9]*","",grpID)
  alnFile <- file.path(codonAlnPath,paste0(grpID.noSplit,".cds.aln"))
  
  # convert sequence alignment to phylip format
  # Only include the sequences from the split tree
  seqIDs.splitClan <- tree$tip.label
  alnPhylip <- readFasta2phylip(alnFile, seqIDs.splitClan)

  writeLines(alnPhylip, "aln.phylip") # temp file
  
  
  # patterns for recognizing the species sets
  outGrpPattern <- "^Os|^Sb|^Zm"
  
  spcPatterns <- c(
    NaSt="^NaSt",
    StLa="^StLa",
    MeNu="^MeNu",
    BrDi="^B",
    HoVu="^H",
    LoPe="^LoPe"
  )
  
  spc2pattern <- function(spcs){
    paste(spcPatterns[spcs],collapse="|")
  }
  
  hSpcs <- list(
    H4c = c("LoPe","HoVu"),
    H4b = c("LoPe","HoVu","BrDi"),
    H4e = c("LoPe","HoVu","BrDi","MeNu"),
    H4d = c("LoPe","HoVu","BrDi","MeNu","StLa"),
    H4a = c("LoPe","HoVu","BrDi","MeNu","StLa","NaSt")
  )

  # first branching species after split exists
  hasFirstBranchingSpc <- function(tree,h){
    firstBranching <- rev(hSpcs[[h]])[1]
    otherSubSpcs <- rev(hSpcs[[h]])[-1]
    return( sum(grepl(spc2pattern(firstBranching),tree$tip.label)) > 0 &
            sum(grepl(spc2pattern(otherSubSpcs),tree$tip.label)) > 0  )
  }

  hasFirstBranchingSplit <- function(tree,h){
    otherSubSpcs <- rev(hSpcs[[h]])[-1]
    return( isClan(tree,grepl(spc2pattern(hSpcs[[h]]),tree$tip.label)) &
              isClan(tree,grepl(spc2pattern(otherSubSpcs),tree$tip.label))  )
  }
  
  # reroot the tree with the outGroup
  tree <- root( phy = tree, outgroup = grep(outGrpPattern,tree$tip.label) )
  
  # remove internal node labels from trees
  tree$node.label <- rep("",tree$Nnode)
  # remove edge lengths
  tree$edge.length <- NULL
  
  # for each hypothesis:
  for( h in names(hSpcs)){
    # Only run if subGrpPattern exists in tree and
    # if current split and next split is good
    if( hasFirstBranchingSpc(tree,h) & hasFirstBranchingSplit(tree,h) ){
      # mark the tree
      markedTree <- markTree(tree, spc2pattern(hSpcs[[h]]))
      
      # Check if H4a mark ends up on the root (occurs if there is only one out-species)
      if( h=="H4a" & grepl("#1;$",write.tree( markedTree )) ){
        cat("resolve root for",grpID,"\n")
        
        # resolve the root        
        rootTree <- root( phy = tree, outgroup = grep(outGrpPattern,tree$tip.label),
                          resolve.root = T )
        # remove labels
        rootTree$node.label <- rep("",tree$Nnode)

        # mark the tree again
        markedTree <- markTree(rootTree, spc2pattern(hSpcs[[h]]))
      }
      
      #plot.phylo(markedTree,show.node.label = T,main=h)
      
      runCodeML( markedTree=markedTree,
                 outFile=file.path( outGrpDir, paste(grpID,h,"H0.out",sep = "_") ),
                 fix_omega = 1 )    
    
      runCodeML( markedTree=markedTree,
                 outFile=file.path( outGrpDir, paste(grpID,h,"H1.out",sep = "_") ),
                 fix_omega = 0 )
    }
  }
}

# main loop:
#
# get the good trees with good topology
# run codeml for each tree inside a temporary folder
mainLoopPAML <- function(x, goodTreesFile, goodTreeStatFile, codonAlnPath, arraySize){
  
  # read trees
  goodTrees <- readRDS(goodTreesFile)
  # read tree metadata
  goodTreeStats <- readRDS(goodTreeStatFile)
  
  # check if core pooids form a clan and each species form separate clans
  hasGoodTopology <- goodTreeStats$isCoreClan & goodTreeStats$allSpcAreClans
  goodTopoTrees <- rev(goodTrees[hasGoodTopology]) # reverse to get smallest first
  
  outDir <- getwd() # output to job's working directory
  
  
  # divide the trees among the jobs in the job array
  for(i in seq(from = x, to = length(goodTopoTrees), by = arraySize)){
    
    tree <- goodTopoTrees[[i]]
    grpID <- names(goodTopoTrees)[i]
  
    # check if Rice is in outGroup
    if( length(grep("^Os",tree$tip.label)) == 0 ){
      cat(grpID,"no rice!?! next!\n")
      next # no rice!?! next!
    }
    
    # check out-group topology
    if( sum(grepl("^Sb|^Zm",tree$tip.label))!=0 & !isClan(tree,!grepl("^Sb|^Zm",tree$tip.label)) ){
      cat(grpID,"Rice don't form a clade with pooids! next!\n")
      next # Rice don't form a clade with pooids! next!
    }

    cat(grpID,"\n")
    
    # create a temp folder
    tmpDir <- tempfile(pattern="codeml")
    dir.create(tmpDir)
    
    # change working directory to temp
    setwd(tmpDir)
    
    # run codeml on tree[i]
    # create a directory for each grp
    outGrpDir <- file.path(outDir,grpID)
    dir.create(outGrpDir,showWarnings = F)
    
  
    doPAML(tree, grpID, codonAlnPath, outGrpDir)
    
    # change working directory back
    setwd(outDir)

    # remove temp files
    unlink(tmpDir,recursive = T)  
    
  }
}

