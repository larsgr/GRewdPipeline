library(ape)
library(phangorn)
#
# splitTreesToGrps(treePath, outGrpFile, outSpcs = c("Os_R", "Sb_R", "Zm_R"))
# ----------------
#
# Looks for "*.tree" files in directory "treePath". Files must have pattern "grp12345.*.tree"
# as the grp ID is used.
#
# Checks trees if they need to be split and splits them into sub-groups so that each 
# new group contains one out/in-species split. Out-species are defined by "outSpcs"
#
# Sub-groups are written to "outGrpFile" in the orthoMCL groups.txt format.
#
splitTreesToGrps <- function(treePath, outGrpFile, outSpcs = c("Os_R", "Sb_R", "Zm_R")){
  
  
  #####
  # Define functions:
  
  # Find all clans in the tree that contain exactly one in-clan and one out-clan
  clanFinder = function(tre, 
                        minSize=3, 
                        ut,
                        output=c("tree", "tips"),
                        summaryFUN=NULL) { 
    
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
  
  
  # Check if tips (given by logical vector isInGrp) form a clan in the tree
  isClan <- function(tree, isInGrp){
    clans = getClans(tree)
    any(apply(clans,1,function(clan){
      all(clan==isInGrp)
    }))
  }
  
  
  
  
  
  #####
  # load trees
  treeFiles <- dir(treePath,pattern="\\.tree$")
  treeFiles <- setNames(file.path(treePath,treeFiles),
                        stringr::str_extract(treeFiles,"grp[0-9]+"))
  
  allTrees <- lapply(treeFiles, read.tree)
  
  #####
  # Select trees to be split
  sapply(allTrees,function(tree){
    spcs <- sapply(strsplit(tree$tip.label, "|", fixed=T), "[", 1)
    uSpcs <- unique(spcs)
    any(uSpcs %in% outSpcs) & # has at least one out-species
      sum(!(uSpcs %in% outSpcs)) >= 2 &  # has at least 2 in-species
      !isClan(tree,spcs %in% outSpcs) # don't have in/out split
  }) -> mustBeSplit
  
  
  #####
  # Split trees with clanFinder
  splitTrees <- lapply(allTrees[mustBeSplit], clanFinder, ut=outSpcs, minSize=3, output="tips")
  
  
  #####
  # rename the split trees with number sequence
  lapply(splitTrees,function(trees){
    setNames(trees,seq_along(trees))
  }) -> renamedTrees
  
  
  #####
  # unlist the split trees. Tree names should now be grp12345.1, grp12345.2 etc..
  treesSplit <- unlist(renamedTrees,recursive = F)
  
  
  
  #####
  # write orthoMCL groups.txt style file. E.g: "grp12345.1: SPC1|gene000 SPC1|gene001 SPC2|gene000"
  writeLines(text = paste0(names(treesSplit),": ",lapply(treesSplit,paste, collapse=" ")),
             con = outGrpFile)
  
}
