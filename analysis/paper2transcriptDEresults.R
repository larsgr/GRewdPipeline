#
# This script was made to extract DE results for the candidate genes in paper 2
#
# List of candidate transcripts: "indata/paper2transcripts/transcripts.csv"
# DE results are loaded from "/mnt/NOBACKUP/mariansc/share/orthos/DESeq/DE.RData"
#
# Since some HoVu transcripts were missing (not in orthoGrp), the read counts
# for these missing transcripts were loaded from the RSEM counts and DESeq was
# re-run for HoVu.

library(DESeq2)
library(BiocParallel)
library(tidyverse)
library(stringr)

load("/mnt/NOBACKUP/mariansc/share/orthos/DESeq/DE.RData")

geneTbl <- read_csv("indata/paper2transcripts/transcripts.csv")

geneTbl %>% 
  mutate(
    baseMean = map2_dbl( .x = species, .y = transcript,
                                  ~ DE[[.x]]$resRamp[.y,"baseMean"]),
    LT_log2FoldChange = map2_dbl( .x = species, .y = transcript,
                                  ~ DE[[.x]]$resRamp[.y,"log2FoldChange"]),
    LT_pvalue = map2_dbl( .x = species, .y = transcript,
                          ~ DE[[.x]]$resRamp[.y,"pvalue"]),
    LT_padj = map2_dbl( .x = species, .y = transcript,
                        ~ DE[[.x]]$resRamp[.y,"padj"]),
    ST_log2FoldChange = map2_dbl( .x = species, .y = transcript , 
                                  ~ DE[[.x]]$resPeak[.y,"log2FoldChange"]),
    ST_pvalue = map2_dbl( .x = species, .y = transcript,
                          ~ DE[[.x]]$resPeak[.y,"pvalue"]),
    ST_padj = map2_dbl( .x = species, .y = transcript,
                        ~ DE[[.x]]$resPeak[.y,"padj"])
  ) -> geneTbl


###
# Run DEseq for HoVu

# load counts for genes in orthoGrps
countData <- read.table("/mnt/NOBACKUP/mariansc/share/orthos/exprTbls/HoVu_expected_countTbl.txt", stringsAsFactors = F)

# fix the T-1 name
names(countData) <- sub("T\\.1","T-1",names(countData))


# get the missing genes
missingGenes <-
  geneTbl %>% 
  filter( species=="HoVu") %>%
  filter( !(transcript %in% rownames(countData))) %>%
  with(transcript)

# load missing read counts
sapply(names(countData),function(sampleID){
  read_tsv( paste0("/mnt/NOBACKUP/mariansc/share/RSEM/HoVu/",sampleID,".genes.results")) %>% 
    filter(sub("\\|","_",gene_id) %in% missingGenes ) %>% 
    with( setNames(expected_count,sub("\\|","_",gene_id))[missingGenes])
}) -> missingCounts

# add the missing data
countData <- rbind(countData,missingCounts)
  

# get colData from sample names
colData <- as.data.frame(str_split_fixed(names(countData),"\\.",4)[ ,1:3])
names(colData) <- c("spcPop","timePoint","mix")
TtoF <- c(`T-1`="ramp0",T0="peak0",T1="peak1",T3="ramp1",T4="ramp1")
colData$Tf <- as.factor(TtoF[as.character(colData$timePoint)])


dds <- DESeqDataSetFromMatrix(countData = round(countData),
                              colData = colData,
                              design = formula(~ Tf ))


dds <- DESeq(dds,parallel = T,BPPARAM = MulticoreParam(workers = 20))

# contrast tests
resRamp <- results(dds, contrast=c("Tf","ramp1","ramp0"))
resPeak <- results(dds, contrast=c("Tf","peak1","peak0"))

# join the results
left_join(
  as.data.frame(resPeak[missingGenes, ]) %>%
    transmute(transcript = missingGenes,
              baseMean = baseMean,
              ST_log2FoldChange = log2FoldChange, 
              ST_pvalue = pvalue,
              ST_padj = padj),
  as.data.frame(resRamp[missingGenes, ]) %>%
    transmute(transcript = missingGenes,
              LT_log2FoldChange = log2FoldChange, 
              LT_pvalue = pvalue,
              LT_padj = padj),
  by="transcript"
) -> missingGeneTbl

# add the missing results to the geneTbl
geneTbl[match(missingGeneTbl$transcript, geneTbl$transcript),colnames(missingGeneTbl)] <-
  missingGeneTbl

write_csv(geneTbl,"indata/paper2transcripts/transcriptsDEresults.csv")

# vsd <- varianceStabilizingTransformation(dds)
# m <- assay(vsd)
# colnames(m) <- names(countData)
# 
# 
# 
# as_data_frame(DE$NaSt$vst[NsIDs,]) %>% 
#   mutate(geneID=NsIDs) %>% 
#   gather(sampleID, vst, -geneID) %>% 
#   mutate(timepointOld=sub(".*\\.(T-?[0-4])\\..*","\\1",sampleID)) %>% 
#   mutate(timepoint= c(`T-1`="W0",T0="D0",T1="D1",T3="W4",T4="W9")[timepointOld]) %>%
#   mutate(effect=ifelse(grepl("D",timepoint),"0short-term","1long-term")) %>%
#   mutate(treatment=ifelse(grepl("0",timepoint),"0before","1after")) %>%
#   ggplot(aes(x=treatment,y=vst)) +
#   geom_point() +
#   facet_grid(geneID~effect)
