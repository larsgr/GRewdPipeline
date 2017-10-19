## Overview of Rmd files

### Analysis used in the article

* __pooid_phylo_script.Rmd__ Figure 1: Species phylogeny with paleoclimate and species habitat climate data.
* __FigureDE.Rmd__ Figure 2: sample clustering tree, venn diagram and trimatrix plot
* __CompareWithGreenup.Rmd__ Figure 3: Comparison with the Greenup genes. Also contain a detailed analysis of the phylogeny of the 22 OGs that were not HCOG
* __GOanalysisBd.Rmd__ Figure 4a: GO analysis
* __TreeWithCandidateGroups.Rmd__ Figure 4b: positive selection tests

The following scripts generates some of the underlying data for the figures:

* __superGeneModel.Rmd__ DE analysis of HCOGs.
* __treeAnalysis6.Rmd__ Defines which trees that are ok for positive selection tests at each split
* __pamlResultsOverview.Rmd__ Summarise the positive selection test results


### Analysis used in paper 2

* __paper2transcriptDEresults.R__: Extract DE results for the candidate genes in paper 2


### compareLoliumArticle.Rmd (deprecated)

Compares the proportion of DE genes with the results from an article that performs a vernalization transcriptome experiment on two variants of lolium. The cluster comparison turned out to be not very useful (too many clusters). The venn diagram shows that there is a significant overlap between the DE genes in lolium and the genes that are DE in all our species, but it is biased by genes with many paralogs and the method of defining a group significant if atleast one of the paralogs are significant. (an interresting note is that there is little overlap between the two variants, which shows that the expression response to cold can vary a lot even when closely related)

### compareModels.Rmd (deprecated)

Clustering of species based of cold response. K-means clustering of genes based on lfc.

### comparePriestArticle.Rmd (not used)

Compare with a study of brachypodium cold response. There is clear correlation between the cold-response gene clusters in the article and brachypodium in our study. Not so clear for the other species and there is little overlap with the greenup genes. 


### CompareVithRMDavidson.Rmd (deprecated.. not finished?)



### complexTrees.Rmd (deprecated..)


### PoolVsReps

Additional DESeq analysis on the supergene expression matrix where triplicate samples for T0 and T1 in BrDi and MeNu are analysed individually as a traditional DE analysis and compared with results when using the pooled samples approached.


DEanalysis.Rmd
DEanalysis2allGroups.Rmd
DEstrategyRobustness.Rmd
exportSuppData.Rmd
ExtractTreesForSpeciesTree.Rmd
GeneListTest.Rmd
GenesOfInterest3.Rmd
genesOfInterrest.Rmd
genesOfInterrest2.Rmd
genesOfInterrest4.Rmd
GOanalysis.Rmd
orthoGroupStats.Rmd
paralogAnalysis.Rmd
pathwayHeatmap.Rmd
pathwayHeatmap2.Rmd



ProfileSimilarityMatrix.Rmd
readCountOverview.Rmd
replicateBiasAnalysis.Rmd
SlidingMatrixPlot.Rmd
StatsTable.Rmd
SupergeneClustering.Rmd
timeModel11grps.Rmd
timeModel11grps2.Rmd
timeSeriesModel.Rmd
treeAnalysis.Rmd
treeAnalysis2.Rmd
treeAnalysis3.Rmd
treeAnalysis4.Rmd
treeAnalysis5.Rmd
ValidateWithGreenup.Rmd
ValidateWithVigeland.Rmd