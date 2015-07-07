library(phangorn)
library(methods)

## setting parameters:
makeTree <- function( alignedFastaFile,
                      outDir,
                      outFilePrefix = sub("\\.aln$","",basename(alignedFastaFile)),
                      model.test  = F,
                      bootstrap = 100,
                      type = "AA" ) { # Type of sequences ("DNA" or "AA")

  
  aa.phylo <- function(alignedFastaFile, model.test = F){
    cat('\nReading in data...\n')
    
    dat = read.phyDat(alignedFastaFile, type = "AA", format = 'fasta')
    dm = dist.ml(dat, model="JTT")
    tree = NJ(dm)
    
    if(model.test){
      cat('made NJ tree...starting model testing..\n')
      (mt <- modelTest(dat, model=c("JTT", "LG", "WAG")))
      cat('Finished model testing..performing ML estimation of gene tree')
      fitStart = mt$Model[which.min(mt$BIC)]
      cat('\n', fitStart)
      if(is.null(fitStart)) {
        cat('WARNING: Error in model test:  Using default JTT+G+I')
        fitNJ = pml(tree, dat, model='JTT', k=4, inv=.2)
        fitStart = 'JTT+G+I'
      }
    }
    
    if(!model.test) {
      cat('\nUsing default JTT+G+I\n')
      fitNJ = pml(tree, dat, model='JTT', k=4, inv=.2)
      fitStart = 'JTT+G+I'
    }
    
    model.desc = unlist(strsplit(fitStart, '\\+'))
    if(length(model.desc)>1 & model.desc[2]=='G') { k = 4 } else { k=1 }
    if(length(model.desc)>1 & model.desc[2]=='I') { inv = 0.2 } else {inv=0} # check if this really works <<----- Model: JTT , gamma.cat = 4 prop. invariable sites= 0 
    model.tag <<- fitStart
    cat('\nModel:',model.desc[1], ', gamma.cat =',  k, 'prop. invariable sites=', inv, '\n\n')
    fitNJ = pml(tree, dat, model=model.desc[1], k=k, inv=inv)
    fit = optim.pml(fitNJ, optNni=TRUE, optInv=k>1, optGamma=inv>0)
    
    return(fit)
  }

  dna.phylo <- function(alignedFastaFile, model.test = F){
    cat('\nReading in data...\n')
    
    dat = read.phyDat(alignedFastaFile, type = "DNA", format = 'fasta')
    dm = dist.ml(dat)
    tree = NJ(dm)
    
    if(model.test){
      cat('Starting model testing..\n')
      mt <- modelTest(dat)
      model.tag <<- mt$Model[which.min(mt$BIC)]
      if(is.null(model.tag)) {
        cat('WARNING: Error in model test!')
        model.test <- F
      } else {
        cat('Using model:',model.tag)
        
        env = attr(mt, "env")
        fitStart = eval(get(model.tag, env), env)
        fit = optim.pml(fitStart, optNni=TRUE, optGamma=TRUE, optInv=TRUE,
                        model=sub("\\+.*$","",model.tag))        
      }
    }
    
    if(!model.test) {
      cat('\nUsing default GTR+G+I\n')
      model.tag <<- 'GTR+G+I'
      fitStart = pml(tree, dat, k=4, inv=.2)
      fit = optim.pml(fitStart, TRUE, TRUE, TRUE, TRUE, TRUE)
    }
    
    return(fit)
  }
  
  
  if( type=="AA"){
    fit <- aa.phylo(alignedFastaFile, model.test = F)
  } else if( type=="DNA"){
    fit <- dna.phylo(alignedFastaFile, model.test = F)    
  } else {
    stop("Unknown type parameter: ",type)
  }
  
  
  options(warn=1)
  
  if( bootstrap > 0 & length(fit$tree$tip.label) < 3){
    cat('\nToo few sequences to bootstrap...\n')
    bootstrap <- 0
  }
  
  if( bootstrap==0 ){
    file.out = file.path(outDir, paste0(outFilePrefix, '_', model.tag, '.tree'))
    cat('\nSaving results to', file.out, '\n')
    write.tree(fit$tree, file=file.out)
  } else {
    cat('\nMaking', bootstrap,  'bootstrap reps....\n')
    fit.bs  <- bootstrap.pml(fit, bs = bootstrap, trees = TRUE,  optNni=TRUE)
    cat('\nBootstrapping done!')
    file.out = file.path(outDir, paste0(outFilePrefix, '_', model.tag, '_BS', bootstrap , '.tree'))
    bs.tree <- plotBS(fit$tree, fit.bs)
  
    cat('\nSaving results to', file.out, '\n')
    write.tree(bs.tree, file=file.out)
  }
}
