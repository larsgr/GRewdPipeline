library(ape)
library(phangorn)
library(stringr)

#####
# define constants:
outSpcs <- c("Os_R", "Sb_R", "Zm_R")


#####
# Define functions:

# Check if tips (given by logical vector) form a clan in the tree
isClan <- function(tree, isInGrp){
  clans = getClans(tree)
  any(apply(clans,1,function(clan){
    all(clan==isInGrp)
  }))
}

# Check if tips (given by vector of species) form a clan in the tree
isClanSpc <- function(tree, clanSpcs){
  isInGrp <- grepl(paste0("^(",paste(clanSpcs,collapse="|"),")"),tree$tip.label)
  isClan(tree,isInGrp)
}

selectTreesToSplit <- function (allTrees) {
  sapply(allTrees,function(tree){
    spcs <- sapply(strsplit(tree$tip.label, "|", fixed=T), "[", 1)
    uSpcs <- unique(spcs)
    any(uSpcs %in% outSpcs) & # has at least one out-species
      sum(!(uSpcs %in% outSpcs)) >= 2 &  # has at least 2 in-species
      !isClan(tree, spcs %in% outSpcs) # don't have in/out split
  }) 
}

loadTrees <- function (treePath) {
  treeFiles <- dir(treePath,pattern="\\.tree$")
  treeFiles <- setNames(file.path(treePath,treeFiles),
                        str_extract(treeFiles,"grp[0-9]+"))
  
  return( lapply(treeFiles, read.tree) )
}

# Find all clans in the tree that contain exactly one in-clan and one out-clan
clanFinder = function(tre, minSize=3, ut, output=c("tree", "tips")){ 
  
  clans = getClans(tre) # phangorn
  
  # get species names (split on pipe)
  specs = sapply(strsplit(colnames(clans), "|", fixed=T), "[", 1)
  
  # convert clans to other format
  clanList_all = apply(clans, 1, function(r) which(as.logical(r)))
  
  # which clans >= minSize and does not contain any outgroup species?
  clanList_in = clanList_all[ vapply(clanList_all, function(cln){
    length(cln) >= minSize  &&  !any(specs[cln] %in% ut)
  }, logical(1))]
  
  # which clans >= 2 and only contain outgroup species?
  clanList_ut = clanList_all[vapply(clanList_all, function(cln) {
    length(cln) >= 2 && all(specs[cln] %in% ut)
  }, logical(1))]
  
  # for each clan
  #   clan is goodClan if: clan > minSize and 
  #   contains exactly one in-clan (size >= minsize) and one out-clan
  goodClan = vapply(clanList_all, function(cln) {
    if(length(cln) < minSize+1) return(F)
    sp = specs[cln]             # associated species
    utgruppe = cln[sp %in% ut]  # out group clan members
    ingruppe = cln[!sp %in% ut] # in group clan members
    if(length(ingruppe) < minSize || length(utgruppe) == 0) # at least 'minsize' in and 1 out
      return(F)
    in_isClan = list(ingruppe) %in% clanList_in  # test if in group members form a clan
    ut_isClan = length(utgruppe) == 1 || list(utgruppe) %in% clanList_ut # same for out members (if more than 1)
    in_isClan && ut_isClan
  }, logical(1))
  
  tip_seq = seq_along(tre$tip.label)
  switch(match.arg(output),
         tree = lapply(clanList_all[goodClan], function(cln) drop.tip(tre, setdiff(tip_seq, cln))),
         tips = lapply(clanList_all[goodClan], function(cln) colnames(clans)[cln]))
}

# A wrapper that also renames the trees
runClanFinder <- function (trees) {
  unnamedTrees <- lapply(trees, clanFinder, ut=outSpcs, minSize=3, output="tree")
  
  #####
  # rename the split trees with number sequence
  lapply(unnamedTrees,function(trees){
    setNames(trees,seq_along(trees))
  }) -> renamedTrees
  
  # return flat list of split trees  
  return( unlist(renamedTrees,recursive = F) )
}

getMinReqTrees <- function(allTrees, splitTrees) {
  
  # requirement 1: must have outGroup as a clan
  isOutGroupClan <- sapply(allTrees, isClanSpc, outSpcs)
  trees <- c(splitTrees, allTrees[isOutGroupClan])
  
  # requirement 2: have at least 1 basal and 1 de-novo core species
  
  sapply(trees, function(tree){
    (sum(grepl("^(HoVu|BrDi)",tree$tip.label)) > 0) &
      (sum(grepl("^(MeNu|NaSt|StLa)",tree$tip.label)) > 0)
  }) -> has1coreAnd1Basal
  
  return( trees[has1coreAnd1Basal] )  
}

getAdditionalStats <- function(trees){
  # species
  spcs <- c("MeNu1","MeNu2","NaSt","StLa","HoVu","BrDi","Hv_R","Bd_R")
  sapply(spcs, function(spc){
    sapply(trees, function(tree){
      sum(grepl(paste0("^",spc),tree$tip.label))
    })
  })-> nSpcs2
  
  nCoreSpc <- rowSums(nSpcs2[,c("HoVu","BrDi","Hv_R","Bd_R")]>0)
  nCoreDeNovoSpc <- rowSums(nSpcs2[,c("HoVu","BrDi")]>0)
  nBasalSpc <- rowSums(nSpcs2[,c("MeNu1","MeNu2","NaSt","StLa")]>0)
  
  has2BasalAnd2Core <- nCoreSpc>1 & nBasalSpc>1
  hasAllDenovo <- nCoreDeNovoSpc==2 & nBasalSpc==4
  hasAllPooids <- nCoreSpc==4 & nBasalSpc==4
  
  # has core split
  sapply(trees,function(tree){
    isClan(tree,grepl("^(H|B)",tree$tip.label))
  }) -> isCoreClan
  
  hasParalogs <- rowSums(nSpcs2>1)>0
  
  sapply(c("H","B","MeNu","NaSt","StLa"), function(spc){
    sapply(trees,function(tree){
      isInGrp <- grepl(paste0("^",spc),tree$tip.label)
      if( sum(isInGrp) < 2){
        return(TRUE) # return TRUE if tip not existing
      } else {
        return(isClan(tree,isInGrp))      
      }
    })
  }) -> isSpcClan
  
  allSpcAreClans <- apply(isSpcClan,1,all)
  
  return( data.frame( allSpcAreClans = allSpcAreClans,
                      isCoreClan = isCoreClan,
                      has2BasalAnd2Core = has2BasalAnd2Core,
                      hasAllDenovo = hasAllDenovo,
                      hasAllPooids = hasAllPooids,
                      hasParalogs = hasParalogs)
  )
}

# write orthoMCL groups.txt style file. E.g: "grp12345.1: SPC1|gene000 SPC1|gene001 SPC2|gene000"  
writeTreesAsGroups <- function(trees,outFile){
  lapply(trees,function(tree){
    paste(tree$tip.label, collapse=" ")
  }) -> tips
  writeLines(text = paste0(names(trees),": ",tips),
             con = outFile)
}


#####
#
# splitTreesToGrps(treePath)
# ----------------
#
# Looks for "*.tree" files in directory "treePath". Files must have pattern "grp12345.*.tree"
# as the grp ID is used.
#
# Checks trees if they need to be split and splits them into sub-groups so that each 
# new group contains one out/in-species split.
#
# 
#
splitTreesToGrps <- function(treePath){
  
  
  #####
  # load trees
  allTrees <- loadTrees(treePath) 
  
  #####
  # Select trees to be split
  mustBeSplit <- selectTreesToSplit(allTrees)
  
  #####
  # Split trees with clanFinder
  splitTrees <- runClanFinder(allTrees[mustBeSplit])
  
  #####
  # Get the trees that fits the minimum requirements
  goodTrees <- getMinReqTrees(allTrees, splitTrees)

  #####
  # Get more stats about each tree (isCore, isAllSpcClans etc..)
  goodTreeStats <- getAdditionalStats(goodTrees)
  
  #####
  # write files
  
  writeTreesAsGroups(splitTrees,"splitGroups.txt")
  writeTreesAsGroups(goodTrees,"goodGroups.txt")

  saveRDS(goodTrees, "goodTrees.rds")
  saveRDS(goodTreeStats, "goodTreeStats.rds")
  
}
