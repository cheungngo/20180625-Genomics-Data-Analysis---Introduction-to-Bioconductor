# Phenotypes Assessment


source("http://www.bioconductor.org/biocLite.R")
biocLite('COPDSexualDimorphism.data')
library(COPDSexualDimorphism.data)
data(lgrc.expr.meta)


# The variable pkyrs in the expr.meta data.frame represents pack years smoked
# Other variables include gender (self-explanatory) and diagmaj (disease status)

names(expr.meta)
# [1] "tissueid"    "sample_name" "newid"       "GENDER"      "age"         "cigever"     "pkyrs"       "diagmaj"    
# [9] "gender"


# What is the number of female participants in this study?:

unique(expr.meta$gender)
# Levels: 1-Male 2-Female

sum(expr.meta$gender=='2-Female')
# [1] 110


# What is the median of the distribution of pack years smoked in this cohort (women and men)?

median(expr.meta$pkyrs)
# [1] 40


# True or False: The distribution of pack-years smoked is well-approximated by a Gaussian (Normal) probability distribution.

qqnorm(expr.meta$pkyrs)
# FALSE


# Which of the following is an aspect of the display that would suggest caution in using the t test 
# in comparing males and females with respect to pack years smoked?

boxplot(pkyrs~gender, data=expr.meta)
# plot 01
# Distributions appear quite asymmetric, with long tails skewed towards high values.


# Variable transformation using boxcox


# transform the pkyrs into a positive var for analysis
expr.meta$pyp1 = expr.meta$pkyrs+1

lm1 = lm(pyp1~gender, data=expr.meta)
library(MASS)
bc1 = boxcox(lm1)
lambda = bc1$x[which(rank(-bc1$y)==1)]
lambda # [1] 0.4646465
# if lambda is 0.5, we use sqrt(pyp1)


boxplot(I(pyp1^lambda)~gender, data=expr.meta)
# plot 02
# the skewness seems disappeared


###


# Chromosomes and SNPs assessment


library(BiocInstaller)
biocLite("BSgenome.Hsapiens.UCSC.hg19")

library(BSgenome.Hsapiens.UCSC.hg19)
BSgenome.Hsapiens.UCSC.hg19

# We can access chromosome 11 like this:
chr11seq <- BSgenome.Hsapiens.UCSC.hg19[["chr11"]]

# Here, for example, is a segment of 25 bases starting  at base 1 million 
subseq(chr11seq,start=10^6,width=25)
# 25-letter "DNAString" instance
seq: GTTTCACCGTGTTAGCCAGGGTGGT


# Read the help file for the fuction countPattern and tell us 
# which of the following sequences is most common on chromosome 11: "ATG", "TGA", "TAA", and "TAG"

?countPattern
testseq = c('ATG','TGA','TAA','TAG')
sapply(testseq, countPattern, subject = chr11seq)
#     ATG     TGA     TAA     TAG 
# 2389002 2561021 2624324 1689356


# Read the help page for the function alphabetFrequency 
# and use it to determine what percent of chromosome 7 is T,C,G, and A.
# What proportion are Cs (including counts of N in the total) (N for the positions not being called)

?alphabetFrequency
chr7seq = BSgenome.Hsapiens.UCSC.hg19[['chr7']]
alphabetFrequency(chr7seq)
#        A        C        G        T        M        R        W        S        Y        K        V        H        D 
# 45997757 31671670 31636979 46047257        0        0        0        0        0        0        0        0        0 
#        B        N        -        +        . 
#        0  3785000        0        0        0 

alphabetFrequency(chr7seq)[2]/sum(alphabetFrequency(chr7seq))
# 0.1990193


# Locations of SNPs in humans

# SNP information is not on the human genome reference sequence. 
# Instead, this information is stored in databases such as dbSNP.

biocLite('SNPlocs.Hsapiens.dbSNP144.GRCh37')
library(SNPlocs.Hsapiens.dbSNP144.GRCh37)

# checking on SNPs
snps144 = SNPlocs.Hsapiens.dbSNP144.GRCh37
s17 = snpsBySeqname(snps144, "17")
head(s17)
# GPos object with 6 positions and 2 metadata columns:
#       seqnames       pos strand |   RefSNP_id alleles_as_ambig
#          <Rle> <integer>  <Rle> | <character>      <character>
#   [1]       17        52      * | rs556541063                M
#   [2]       17        56      * | rs145615430                Y
#   [3]       17        78      * | rs148170422                S
#   [4]       17        80      * | rs183779916                R
#   [5]       17        92      * | rs562410061                K
#   [6]       17       168      * | rs529798787                R
#   -------
#   seqinfo: 25 sequences (1 circular) from GRCh37.p13 genom

# The first one listed is rs556541063 which is at location 52.
# What is the location on chr17 of SNP rs73971683?

snpsById(snps144,'rs73971683')
# 135246


# GWAS: Linking SNP genotypes to disease risk

# The Bioconductor gwascat package includes information on a catalog of GWAS results assembled at EMBL-EBI

biocLite('gwascat')
library(gwascat)
data(ebicat37)
ebicat37
# gwasloc instance with 22688 records and 36 attributes per record.
# Extracted:   
# Genome:  GRCh37 
# Excerpt:
# GRanges object with 5 ranges and 3 metadata columns:
#    seqnames                 ranges strand |                  DISEASE/TRAIT        SNPS   P-VALUE
#       <Rle>              <IRanges>  <Rle> |                    <character> <character> <numeric>
# [1]    chr11 [ 41820450,  41820450]      * | Post-traumatic stress disorder  rs10768747     5e-06
# [2]    chr15 [ 35060463,  35060463]      * | Post-traumatic stress disorder  rs12232346     2e-06
# [3]     chr8 [ 97512977,  97512977]      * | Post-traumatic stress disorder   rs2437772     6e-06
# [4]     chr9 [100983826, 100983826]      * | Post-traumatic stress disorder   rs7866350     1e-06
# [5]    chr15 [ 54715642,  54715642]      * | Post-traumatic stress disorder  rs73419609     6e-06

# check out the @ elementMetadat @ listdata

# Which chromosome has the most GWAS hits in the catalog? Use an integer

# counting the no of chr in gwas
sort(table(ebicat37$CHR_ID),decreasing=TRUE)

# or
which(rank(-tab_gwas)==1)
# 6  # chr6
# 20 # position 20


# You can use the notation mcols(ebicat37)[,"DISEASE/TRAIT"] 
# to get a vector of names of diseases with genetic associations recorded in the gwascat

?mcols
# mcols(x, use.names=FALSE), mcols(x) <- value: Get or set the metadata columns. 
# If use.names=TRUE and the metadata columns are not NULL, 
# then the names of x are propagated as the row names of the returned DataTable object.

# What is the disease/trait with the most associations?

head(sort(table(ebicat37$'DISEASE/TRAIT'),decreasing=T))
# Obesity-related traits                 Height      IgG glycosylation        Type 2 diabetes   Rheumatoid arthritis 
#                    957                    822                    699                    340                    294 
#        Crohn's disease 
#                    249 


###


# Gene Expression Assessment

library(BiocInstaller)
biocLite("genomicsclass/tissuesGeneExpression")
library(tissuesGeneExpression)

data(tissuesGeneExpression)
head(e[,1:5])
table(tissue)


# you have object e, tab, and tissue in your workspace.  
# You can work with them separately 
# but it is preferable in Bioconductor to unify them in an object.  
# In 2017, the preferred object type is SummarizedExperiment.  

biocLite('SummarizedExperiment')
library(SummarizedExperiment)
tissSE = SummarizedExperiment(list(rma=e))
colData(tissSE) = DataFrame(tab) 
# note that 'tab' was arranged by sampleID by rows
# And other qualities of tab would also be arranged into the tissSE


# Localization of expression to tissues

dim(assay(tissSE))
# [1] 22215   189  # genes by row; samples by column
# using the assay() command
# note that the assay() command is also inside the SummarizedExperiment package

# Locating features with ID "209169_at"
mean(assay(tissSE["209169_at",]))
# [1] 7.26365

# stratify assay data for the feature by tissue and create boxplots
boxplot(assay(tissSE)["209169_at",]~tissSE$Tissue, main='plot 03')
# plot 03


# Which of the following ID(s) appear to represent a gene specific to placenta

IDs = c("201884_at", "209169_at", "206269_at", "207437_at", "219832_s_at", "212827_at")

par(mfrow=c(2,3))
sapply(IDs,function(i) {
  boxplot(assay(tissSE)[i,]~tissSE$Tissue)
})
# plot 04


# Discovery of microarray annotation in Bioconductor

# How would we go about finding more information about gene "206269_at" for example? 
# Does it have a known function? Where is it on the genome? What is its sequence?

# The microarray product used to make the measurements described here is the Affymetirx Human GeneChip HG133A. 
# Search the Bioconductor website and determine which of the following packages 
# provides a connection to gene information:

# hgu133a.db 

biocLite('hgu133a.db')


# Oligo sequences on affymetrix arrays

# The affymetrix technology for mRNA abundance measurement is based on 
# hybridization of highly processed biological sample material to synthetic oligonucleotide "probes" 
# that are on fixed positions of the microarray surface. 
# Bioconductor provides detailed information on the probe and array structure as published by affymetrix.

# Install and attach the hgu133aprobe package.

biocLite('hgu133aprobe')
library(hgu133aprobe)
head(hgu133aprobe)
#                    sequence   x   y Probe.Set.Name Probe.Interrogation.Position Target.Strandedness
# 1 CACCCAGCTGGTCCTGTGGATGGGA 467 181      1007_s_at                         3330           Antisense
# 2 GCCCCACTGGACAACACTGATTCCT 531 299      1007_s_at                         3443           Antisense
# 3 TGGACCCCACTGGCTGAGAATCTGG  86 557      1007_s_at                         3512           Antisense
# 4 AAATGTTTCCTTGTGCCTGCTCCTG 365 115      1007_s_at                         3563           Antisense
# 5 TCCTTGTGCCTGCTCCTGTACTTGT 207 605      1007_s_at                         3570           Antisense
# 6 TGCCTGCTCCTGTACTTGTCCTCAG 593 599      1007_s_at                         3576           Antisense
dim(hgu133aprobe)
# [1] 247965      6

# finding the information of '206269_at'
ind = which(hgu133aprobe$Probe.Set.Name=='206269_at')
hgu133aprobe[ind,]
#                        sequence   x   y Probe.Set.Name Probe.Interrogation.Position Target.Strandedness
# 63835 GCTAACCTATCTAGGACCTGATCTA 405 315      206269_at                         2080           Antisense
# 63836 ACCTGATCTATGCCTTCTTGGGAAC 201  91      206269_at                         2095           Antisense
# 63837 AAAGCTGAAGATTCTACCACTGAAG  63 123      206269_at                         2170           Antisense
# 63838 TGACAAGGAACCTGACCTAACCCCA 547 321      206269_at                         2269           Antisense
# 63839 TGACCTAACCCCATTTTTCATAAAG 184  87      206269_at                         2281           Antisense
# 63840 TTCCCCTTGCACAGGGGCTTCTGTT 168 697      206269_at                         2402           Antisense
# 63841 ATCCCATTTATTATTTGTGGCACCT  21  29      206269_at                         2429           Antisense
# 63842 GTGGCACCTATATCAATGTGGGGTT 410 485      206269_at                         2445           Antisense
# 63843 GAAGCTCTTGTGCAATTAGCCTGGC 161 355      206269_at                         2498           Antisense
# 63844 TTAGCCTGGCAATGATTGGCCGGCC 333 687      206269_at                         2513           Antisense
# 63845 GGCCGGCCATGTTTGAATGCCAGTT 549 547      206269_at                         2530           Antisense

# How many oligos are used to interrogate samples for gene GCM1, annotated to probe 206269_at?
length(ind)
# [1] 11


# Annotation enhancement of a SummarizedExperiment

library(hgu133a.db)
sym = mapIds(hgu133a.db, keys=rownames(tissSE), column="SYMBOL", keytype="PROBEID")
nm = mapIds(hgu133a.db, keys=rownames(tissSE), column="GENENAME", keytype="PROBEID")
rowData(tissSE) = DataFrame(symbol=sym, genename=nm)
# note the changes in @ elementMetadat @ listData
# note that $ genename allows us to search with the gene name

# access the gene with genename
tissSE[grep('kinase', rowData(tissSE)$genename),]
# class: SummarizedExperiment 
# dim: 1064 189 
# ...

# How many features are annotated to genes with 'kinase' in their name?
# 1064


###


# Microarray Assessment


# Using a SummarizedExperiment


# using hgu133a again for the tissuesGeneExpression
library(tissuesGeneExpression)
data(tissuesGeneExpression)
library(SummarizedExperiment)
tissSE = SummarizedExperiment(list(rma=e))
colData(tissSE) = DataFrame(tab)
library(hgu133a.db)
sym = mapIds(hgu133a.db, keys=rownames(tissSE), column="SYMBOL", keytype="PROBEID")
nm = mapIds(hgu133a.db, keys=rownames(tissSE), column="GENENAME", keytype="PROBEID")
rowData(tissSE) = DataFrame(symbol=sym, genename=nm)

# How many features in this SummarizedExperiment measure expression of gene H2AFX?
# Note that H2AFX is a symbol

tissSE[grep('H2AFX', rowData(tissSE)$symbol),]
# rownames(4): 205436_s_at 212524_x_at 212525_s_at 213344_s_at

# another code for that
sum(rowData(tissSE)$symbol=="H2AFX", na.rm=TRUE)


# Verify that 205436_s_at is the affymetrix code for H2AFX and then consider the following plot:

par(las=2, mar=c(10,4,2,2))
boxplot(as.numeric(assay(tissSE["205436_s_at",]))~tissSE$Tissue)
# plot 05

# Which of the following relationships are suggested by this plot?
# Expression of H2AFX is greater in hippocampus than in cerebellum. 


###

# Section 2: Structure and management of genome-scale data with Bioconductor  

# Structures tailored to microarray experiments

###


# ExpressionSet


library(GSE5859Subset)
data(GSE5859Subset)
dim(geneExpression) # [1] 8793   24
dim(geneAnnotation) # [1] 8793    4
dim(sampleInfo) # [1] 24  4


# checking if it is all equal

all.equal(rownames(geneExpression),geneAnnotation$PROBEID) # [1] TRUE
all.equal(colnames(geneExpression),sampleInfo$filename) # [1] TRUE


# building an expression set

rownames(geneAnnotation) = geneAnnotation$PROBEID
rownames(sampleInfo) = sampleInfo$filename


# using biobase

library(Biobase)
es5859 = ExpressionSet(assayData = geneExpression)
pData(es5859) = geneAnnotation
fData(es5859) = sampleInfo
es5859
# ExpressionSet (storageMode: lockedEnvironment)
# assayData: 8793 features, 24 samples 
#   element names: exprs 
# protocolData: none
# phenoData
#   sampleNames: 1007_s_at 1053_at ... AFFX-r2-P1-cre-5_at (8793 total)
#   varLabels: PROBEID CHR CHRLOC SYMBOL
#   varMetadata: labelDescription
# featureData
#   featureNames: GSM136508.CEL.gz GSM136530.CEL.gz ... GSM136572.CEL.gz (24 total)
#   fvarLabels: ethnicity date filename group
#   fvarMetadata: labelDescription
# experimentData: use 'experimentData(object)'
# annotation:


# methods for the class 'ExpressionSet'

methods(class = ExpressionSet)
#  [1] $                $<-              [                [[               [[<-             abstract        
#  [7] annotation       annotation<-     as.data.frame    as.matrix        assayData        assayData<-     
# [13] classVersion     classVersion<-   coerce           combine          description      description<-   
# [19] dim              dimnames         dimnames<-       dims             esApply          experimentData  
# [25] experimentData<- exprs            exprs<-          fData            fData<-          featureData     
# [31] featureData<-    featureNames     featureNames<-   fvarLabels       fvarLabels<-     fvarMetadata    
# [37] fvarMetadata<-   initialize       isCurrent        isVersioned      makeDataPackage  notes           
# [43] notes<-          pData            pData<-          phenoData        phenoData<-      preproc         
# [49] preproc<-        protocolData     protocolData<-   pubMedIds        pubMedIds<-      rowMedians      
# [55] rowQ             sampleNames      sampleNames<-    show             storageMode      storageMode<-   
# [61] updateObject     updateObjectTo   varLabels        varLabels<-      varMetadata      varMetadata<-   
# [67] write.exprs     

# The most important methods are
# exprs(): get the numerical expression values
# pData(): get the sample-level data
# fData(): get feature-level data
# annotation(): get a tag that identifies nomenclature for feature names
# experimentData(): get a MIAME-compliant metadata structure


es5859[1:4,1:3]

# ...
# assayData: 4 features, 3 samples 
# ...
# phenoData
#   sampleNames: 1007_s_at 1053_at 117_at
# ...
# featureData
#   featureNames: GSM136508.CEL.gz GSM136530.CEL.gz GSM136517.CEL.gz GSM136576.CEL.gz
# ...


# restricting the ExpressionSet to e.g. chrY

names(pData(es5859))
# [1] "PROBEID" "CHR"     "CHRLOC"  "SYMBOL"

ind = which(pData(es5859)$CHR == 'chrY')
es5859[ind,]
# ExpressionSet (storageMode: lockedEnvironment)
# assayData: 21 features, 24 samples
# ...


biocLite('GEOquery')
library(GEOquery)
GEO5859 = getGEO('GSE5859')
GEO5859
# we got their data in the format of ExpressionSet
# this time 208 samples instead of 24


# experimentData

biocLite('annotate')
library(annotate)
mi = pmid2MIAME("17206142")
experimentData(es5859) = mi
experimentData(es5859)
# Experiment data
#   Experimenter name: Spielman RS 
#   Laboratory: NA 
#   Contact information:  
#   Title: Common genetic variants account for differences in gene expression among ethnic groups. 
#   URL:  
#   PMIDs: 17206142 
#   Abstract: A 145 word abstract is available. Use 'abstract' method.

abstract(experimentData(es5859))
# [1] "Variation in DNA sequence contributes to individual differences in quantitative traits,
# ...


###


# The GEOquery package: ExpressionSets from NCBI's repository


biocLite('GEOquery')
library(GEOquery)
# GEO for gene expression omnibus
# https://www.ncbi.nlm.nih.gov/geo/


# We have an especial interest in the genomics of glioblastoma and have identified a paper (PMID 27746144) 
# addressing a metabolic pathway whose manipulation may enhance treatment development strategies. 
# Affymetrix Primeview arrays were used, with quantifications available in GEO

# from some searching we got the GSE number: GSE78703
GEO78703 = getGEO('GSE78703')[[1]]  # the [[1]] is someehow very important
GEO78703
# ExpressionSet (storageMode: lockedEnvironment)
# assayData: 49395 features, 12 samples 
# ...


# Digging dip on the data

names(pData(GEO78703))
#  [1] "title"                   "geo_accession"           "status"                  "submission_date"        
#  [5] "last_update_date"        "type"                    "channel_count"           "source_name_ch1"        
#  [9] "organism_ch1"            "characteristics_ch1"     "characteristics_ch1.1"   "molecule_ch1"           
#  ...

pData(GEO78703)$title
#  [1] NHA DMSO 1       NHA DMSO 2       NHA DMSO 3       U87vIII DMSO 1   U87vIII DMSO 2   U87vIII DMSO 3  
#  [7] NHA LXR623 1     NHA LXR623 2     NHA LXR623 3     U87vIII LXR623 1 U87vIII LXR623 2 U87vIII LXR623 3
# 12 Levels: NHA DMSO 1 NHA DMSO 2 NHA DMSO 3 NHA LXR623 1 NHA LXR623 2 NHA LXR623 3 U87vIII DMSO 1 ... U87vIII LXR623 3


# Creating a table for some of the phenoData

table(pData(GEO78703)$'cell type:ch1', pData(GEO78703)$'treated with:ch1')
#                                  DMSO control LXR-623 5 uM for 24 hr
#   Normal human astrocytes                   3                      3
#   U87EGFRvIII glioblastoma cells            3                      3


###


# Navigating EMBL's ArrayExpress repository

# https://www.ebi.ac.uk/arrayexpress/


biocLite('ArrayExpress')
library(ArrayExpress)
sets = queryAE(keywords = "glioblastoma", species = "homo+sapiens")
dim(sets)
# [1] 535   8
# same as what we got from the web


head(sets[,-c(7,8)])
#                      ID Raw Processed ReleaseDate PubmedID                     Species
# E-MTAB-4004 E-MTAB-4004 yes        no  2018-06-28 29944140                Homo sapiens
# E-MTAB-6408 E-MTAB-6408  no        no  2018-06-26       NA                Homo sapiens
# E-MTAB-6003 E-MTAB-6003 yes        no  2018-06-12       NA                Homo sapiens
# E-MTAB-6682 E-MTAB-6682  no        no  2018-05-31 29844126 Homo sapiens | Mus musculus
# E-MTAB-6681 E-MTAB-6681  no        no  2018-05-31 29844126                Homo sapiens
# E-MTAB-5552 E-MTAB-5552 yes        no  2018-05-11       NA                Homo sapiens


# then one can access the data by the getAE() command


###


# ExpressionSet assessment


# genefu again, for the breast expression data

library(Biobase)
library(genefu)
data(nkis)
dim(demo.nkis)
head(demo.nkis)[,1:8]
#         dataset series  id er grade node size age
# NKI_123     NKI   NKI2 123  1     3    0  2.0  48
# NKI_327     NKI   NKI2 327  1     2    1  2.0  49
# NKI_291     NKI   NKI2 291  1     2    1  1.2  39
# NKI_370     NKI   NKI2 370  1     1    1  1.8  51
# NKI_178     NKI   NKI2 178  1     2    1  3.0  48
# NKI_176     NKI   NKI2 176  1     2    1  5.0  46

# NKI seems like the samples

head(annot.nkis)

# NM_ are genes


# Try

nkes = ExpressionSet(data.nkis, phenoData=AnnotatedDataFrame(demo.nkis),
                     featureData=AnnotatedDataFrame(annot.nkis))
# 6 errors observed


# After exploring the data, we saw that for data.nkis, rows are samples and cols are genes

nkes = ExpressionSet(assayData=t(data.nkis)) # need to transpose the data.nkis
pData(nkes) = demo.nkis
fData(nkes) = annot.nkis
nkes


# Acquire the dataset associated with the paper 
# "Transformation from committed progenitor to leukaemia stem cell initiated by MLL-AF9", PMID 16862118, 
# by Krivtsov and colleagues. 
# This paper uses Affymetrix microarrays to study how macrophage precursors may become cancer cells.

library(GEOquery)
lstem = getGEO("GSE3725")

# What is the class of lstem after this operation completes?

class(lstem)
[1] "list"

# Transforming back to an ExpressionSet
# How many features and samples are in the ExpressionSet

lstem = lstem[[1]]
lstem
# ExpressionSet (storageMode: lockedEnvironment)
# assayData: 22690 features, 28 samples 

# One common difficulty of working with GEO is that the characteristics of different samples 
# are not always easily determined. 
# Sometimes there is no annotation, and sometimes the annotation is present in an unusual field. 
# In this case, the sample characteristic of interest is the type of cell on which expression measures were taken.

pData(lstem)$title 
# ...
# 28 Levels: GMP expressing GFPrep1 GMP expressing GFPrep2 GMP expressing GFPrep3 ... Lin-, sca+, kit+ (HSC enriched)rep5

# ignoring the first 6
lstem = lstem[, -c(1:6)]

# How many samples are of type L-GMP?
length(grep('L-GMP', pData(lstem)$title))
# [1] 6

# Improving sample and feature labeling for a heatmap

## perform an elementary normalization
ee = exprs(lstem)
ee[ee<0] = 0 
eee = log(ee+1)
## boxplot(data.frame(eee))
meds = apply(eee,2,median)
tt = t(t(eee)-meds)
## boxplot(data.frame(tt))
## assign the normalized values to ExpressionSet
exprs(lstem) = tt

# simplify downstream labeling with gene symbol
featureNames(lstem) = make.names(fData(lstem)$"Gene Symbol", unique=TRUE)
# note the changes in featureNames under featureData

# The following code is somewhat complex, but it simplifies labeling of cell types 
# by stripping away details of marker configurations.

# reformat the naming of cell types
ct = pData(lstem)[,1]
ct = as.character(ct)
cct = gsub(".*(\\(.*\\)).*", "\\1", ct) 
cct = make.unique(cct)
cct = gsub(" enriched", "", cct)
# use the cell types as sample names
sampleNames(lstem) = cct

# Four genes identified in the stemness signature are given in a vector below. 
# We will use these for a small-scale heatmap.
# select some members of the stem cell signature
inds = which(fData(lstem)$"Gene Symbol" %in% c("Stat1", "Col4a1", "Hoxa9", "Itgb5"))

# Obtain a simple heatmap
heatmap(exprs(lstem[inds,]), Colv=NA)

# What's the total number of probes interrogating the four genes of interest?

length(inds)
# [1] 11


###


# IRanges


library(BiocInstaller)
biocLite('IRanges')
library(IRanges)


# start, end, width

ir1 = IRanges(5,10)
ir1
# IRanges object with 1 range and 0 metadata columns:
#          start       end     width
#      <integer> <integer> <integer>
#  [1]         5        10         6

# another way 

ir2 = IRanges(5,width = 6)
ir1 == ir2 # [1] TRUE

# specific commands

start(ir1) # [1] 5
end(ir1) # [1] 10
width(ir1) # [1] 6


# more than a range at a time

ir3 = IRanges(c(5,7),c(8,10))
ir3
#           start       end     width
#       <integer> <integer> <integer>
#   [1]         5         8         4
#   [2]         7        10         4


# intrarange operations

ir4 = shift(ir3,c(0,2))
ir4
#         start       end     width
#     <integer> <integer> <integer>
# [1]         5         8         4
# [2]         9        12         4

ir5 = narrow(ir4[1], start = 2)
ir5
#           start       end     width
#       <integer> <integer> <integer>
#   [1]         6         8         3

# '6' is the 2nd relative to '5'

ir6 = narrow(ir5, end = 2)
ir6
#           start       end     width
#       <integer> <integer> <integer>
#   [1]         6         7         2

# for 'end' in narrow it walks towards the left

ir7 = IRanges(3,20)
ir7
ir8 = flank(ir7,width = 3,start = T,both = F)
ir8
# [1]         0         2         3
ir9 = flank(ir7,width = 3,start = F,both = F)
ir9
# [1]        21        23         3
ir10 = flank(ir7,width = 3,start = T,both = T)
ir10
# [1]         0         5         6
ir11 = flank(ir7,width = 3,start = F,both = T)
ir11
# [1]        18        23         6

# both = T meaning the inclusion of both inside and outside the range

ir12 = IRanges(c(3,5,17),c(10,8,20))
ir12
#   [1]         3        10         8
#   [2]         5         8         4
#   [3]        17        20         4

range(ir12) 
#   [1]         3        20        18
reduce(ir12)
#   [1]         3        10         8
#   [2]        17        20         4
gaps(ir12)
#   [1]        11        16         6
disjoin(ir12)
#   [1]         3         4         2
#   [2]         5         8         4
#   [3]         9        10         2
#   [4]        17        20         4


###


# IRanges assessment


# zooming in with '*'

# Load the IRanges package. Define an integer range starting at 101 and ending at 200. 
# If we use the operation *2, this will zoom in, giving us a range with half the width. 
# What is the starting point of the resulting range?

ir13 = IRanges(101,200)
ir13*2
# [1]       126       175        50


# If we use the operation narrow(x, start=20), what is the new starting point of the range?

narrow(ir13, start = 20)
# [1]       120       200        81


# zooming out with '+'
# If we use the operation +25, what is the width of the resulting range?

ir13+25
# [1]        76       225       150


# Define an IRanges with starts at 1,11,21 and ends at 3,15,27. 
# width() gives the widths for each range. What is the sum of the widths of all the ranges?

ir14 = IRanges(c(1,11,21),c(3,15,27))
sum(width(ir14))
# [1] 15


# Visualizing and projecting ranges

# Define an IRanges object, x, with the following set of ranges:
  
# Starts at 101,106,201,211,221,301,306,311,351,361,401,411,501
# Ends at 150,160,210,270,225,310,310,330,390,380,415,470,510
ir15 = IRanges(c(101,106,201,211,221,301,306,311,351,361,401,411,501),
               c(150,160,210,270,225,310,310,330,390,380,415,470,510))

# Plot these ranges using the plotRanges function in the ph525x package. 
# You can install this library if you have not done so already, with the command: 
biocLite("genomicsclass/ph525x")
library(ph525x)

# What is the total width from 101 to 510 which is not covered by ranges in x?

plotRanges(ir15)
# plot 07

sum(width(gaps(ir15)))
# [1] 130


# How many disjoint ranges are contained within the ranges in 'x' from the previous question? 
# By disjoint ranges, we mean the following: for two ranges [1,10] and [6,15], 
# there are three disjoint ranges contained within: [1,5], [6,10], and [11,15].

disjoin(ir15)
# 17


# resize()

resize(ir15,1)
#    [1]       101       101         1
#    [2]       106       106         1
#    [3]       201       201         1
#    [4]       211       211         1
#    [5]       221       221         1
#    ...       ...       ...       ...
#    [9]       351       351         1
#   [10]       361       361         1
#   [11]       401       401         1
#   [12]       411       411         1
#   [13]       501       501         1


###


# Genomic ranges: GRanges


library(BiocInstaller)
biocLite("GenomicRanges")
library(GenomicRanges)


# Specifying GenomicRanges  # as an extension from IRanges

ir12 = IRanges(c(3,5,17),c(10,8,20))
gr12 = GRanges(seqnames = 'chrZ', ranges = ir12, strand = '+', seqlengths = c(chrZ=100))
gr12
# GRanges object with 3 ranges and 0 metadata columns:
#       seqnames    ranges strand
#          <Rle> <IRanges>  <Rle>
#   [1]     chrZ  [ 3, 10]      +
#   [2]     chrZ  [ 5,  8]      +
#   [3]     chrZ  [17, 20]      +
#   -------
#   seqinfo: 1 sequence from an unspecified genome


# Operating on GRanges

shift(gr12,2)
#   [1]     chrZ  [ 5, 12]      +
#   [2]     chrZ  [ 7, 10]      +
#   [3]     chrZ  [19, 22]      +

shift(gr12,90)
# error (out of bound)

trim(shift(gr12,90))
#   [1]     chrZ [ 93, 100]      +
#   [2]     chrZ [ 95,  98]      +
#   [3]     chrZ [101, 100]      +

# Other operations are the same as the IRanges
# e.g. narrow(), flank(), range(), reduce(), gaps(), disjoin()


# metadata columns
# adding columns to the Granges Object

mcols(gr12)$values = c(1,2,3)
gr12
# GRanges object with 3 ranges and 1 metadata column:
#       seqnames    ranges strand |    values
#          <Rle> <IRanges>  <Rle> | <numeric>
#   [1]     chrZ  [ 3, 10]      + |         1
#   [2]     chrZ  [ 5,  8]      + |         2
#   [3]     chrZ  [17, 20]      + |         3


# GrangesList

ir12 = IRanges(c(3,5,17),c(10,8,20))
gr12 = GRanges(seqnames = 'chrZ', ranges = ir12)
ir14 = IRanges(c(1,11,21),c(3,15,27))
gr14 = GRanges(seqnames = 'chrZ', ranges = ir14)
gr_12_14 = GRangesList(gr12,gr14)
gr_12_14
# GRangesList object of length 2:
# [[1]] 
# GRanges object with 3 ranges and 0 metadata columns:
# seqnames    ranges strand
# <Rle> <IRanges>  <Rle>
#   [1]     chrZ  [ 3, 10]      *
#   [2]     chrZ  [ 5,  8]      *
#   [3]     chrZ  [17, 20]      *
# [[2]] 
# GRanges object with 3 ranges and 0 metadata columns:
# seqnames   ranges strand
# [1]     chrZ [ 1,  3]      *
# [2]     chrZ [11, 15]      *
# [3]     chrZ [21, 27]      *
# -------
# seqinfo: 1 sequence from an unspecified genome; no seqlengths

length(gr_12_14)
# [1] 2

gr_12_14[[1]]
# allow access of the first GRanges object


# Finding overlaps in GRanges objects

fo_12_14 = findOverlaps(gr12,gr14)
fo_12_14
# Hits object with 1 hit and 0 metadata columns:
#     queryHits subjectHits
#     <integer>   <integer>
# [1]         1           1
# -------
# queryLength: 3 / subjectLength: 3

# first element of gr12 overlaps with the first element of gr14

queryHits(fo_12_14) # [1] 1
subjectHits(fo_12_14) # [1] 1

gr12 %over% gr14
# [1]  TRUE FALSE FALSE
# note that this is with respect to the first object i.e. gr12

countOverlaps(gr12,gr14)
# [1] 1 0 0


# Rle
# Run length encoding

r = Rle(c(-1,-1,1,1,0,0,0,rep(-1,20))) # note the capital R in 'Rle'
r
# numeric-Rle of length 27 with 4 runs
# Lengths:  2  2  3 20
# Values : -1  1  0 -1

Views(r,c(5,7),c(12,15))
# Views on a 27-length Rle subject
# views:
#     start end width
# [1]     5  12     8 [ 0  0  0 -1 -1 -1 -1 -1]
# [2]     7  15     9 [ 0 -1 -1 -1 -1 -1 -1 -1 -1]


###


# GRanges Assessment


# Understanding strand orientation with resize

# In the first week, in the subsection "What We Measure and Why", 
# we learned that DNA has two strands. 
# These two strands are often called plus, "+", and minus, "-".

# The GRanges object in the GenomicRanges package extends the concept of interval ranges 
# in two major ways. The ranges are now also identified by:
# 1. the chromosome we are referring to 
# (in Bioconductor, this is called "seqnames")
# 2. the strand of the DNA we are referring to ("+" or "-"). 
# No strand is labelled with a star, "*".

# Without these two pieces of information, 
# a specification of a range of DNA would be ambiguous. 
# Let's make two ranges, with strand and chromosome information, 
# and see how the range operations act based on strand.

x = GRanges("chr1", IRanges(c(1,101),c(50,150)), strand=c("+","-"))

# In the last assessment, we visualized IRanges 
# with the plotRanges function in the ph525x library. 
# We can get the internal IRanges from a GRanges object with the following code:

ranges(x)
#   [1]         1        50        50
#   [2]       101       150        50

# So let's define a new plotting function:
  
library(ph525x)
plotGRanges = function(x) plotRanges(ranges(x))
# one cannot simply plot GRanges object using ranges()

# Compare x and resize(x,1) using plotGRanges. 
# The result of running resize(x,1) is two ranges of width 1 which start...

x_1 = resize(x,1)
par(mfrow=c(2,1))
plotGRanges(x);plotGRanges(x_1)
# plot 08

# at the left-most point of the "+" strand ranges in x, 
# and the right-most point of the "-" strand ranges in x


# Intersecting transcripts with basic operations

# Suppose we have two different sets of ranges, 
# which overlap somewhat but not entirely. 
# This is the case for many genes, 
# in which there are different versions of transcripts, also called isoforms. 
# The different transcripts consist of exons which end up 
# in the final mRNA molecule, and a number of transcripts 
# can share exons or have exons which are overlapping but not identical ranges.

# We'll start with a toy example, and learn how to load real genes later:

x = GRanges("chr1", IRanges(c(101,201,401,501),c(150,250,450,550)), strand="+")

y = GRanges("chr1", IRanges(c(101,221,301,401,541),c(150,250,350,470,550)), 
            strand="+")

# Plot these two sets of ranges using par(mfrow=c(2,1)) 
# and two calls to plotGRanges.

par(mfrow=c(2,1))
plotGRanges(x);plotGRanges(y)
# plot 09

# If we want to keep the information about which set the ranges belong to, 
# we could combine the two GRanges into a GRangesList:

GRangesList(x,y)

# However, if we want to combine them into a single GRanges, we can use c():

c(x,y)

# Find the total width which is covered by ranges in both x and y. 
# Hint: use c(), disjoin() and %over%.

g_xy = c(x,y)
g_xy
d_xy = disjoin(g_xy)
ind = d_xy %over% x & d_xy %over% y
o_xy = ranges(d_xy[ind])
sum(width(o_xy)) # [1] 140


# Subregions that distinguish transcripts

# What is the total width which is in x or y but not in both?

sum(width(ranges(d_xy[!ind])))
# [1] 130


# The role of strand labeling in range operations

# Define a new genomic range, 'z', which covers range(ranges(x)) 
# but has the opposite strand.

z = GRanges('chr1',ranges(x),strand='-')

# What is the number of ranges in x which overlap z 
# according to the %over% command?

sum(x %over% z)
# [1] 0


###


# Operating on GRanges

library(IRanges)
ir <- IRanges(c(3, 8, 14, 15, 19, 34, 40),
              width = c(12, 6, 6, 15, 6, 2, 7))

library(GenomicRanges)


# amending data in GRanges object

gir = GRanges('chr1',ir)
strand(gir) = c(rep('+',4),rep('-',3))
gir
#       seqnames    ranges strand
#          <Rle> <IRanges>  <Rle>
#   [1]     chr1  [ 3, 14]      +
#   [2]     chr1  [ 8, 13]      +
#   [3]     chr1  [14, 19]      +
#   [4]     chr1  [15, 29]      +
#   [5]     chr1  [19, 24]      -
#   [6]     chr1  [34, 35]      -
#   [7]     chr1  [40, 46]      -


# plotting GRanges

par(mfrow=c(4,1))

# original

library(ph525x)
plotRanges(ranges(gir),xlim=c(0,60),main='original')

# the starting point of the gene

plotRanges(ranges(resize(gir,1)),xlim=c(0,60),main='starting point')
# note that the starting point would be different for '+' and '-' strands

# the promoter (supposedly 3 basepairs before each gene)

plotRanges(ranges(flank(gir,3,start=T)),xlim=c(0,60),main='promotor')

# downstream promotor with the length of two

plotRanges(ranges(flank(gir,2,start = F)),xlim = c(0,60),main='downstream promotor')

# plot 10


###


# Finding Overlaps


library(BiocInstaller)
biocLite('GenomicFeatures')
biocLite('genomicsclass/ERBS')

library(GenomicFeatures)
library(ERBS)

browseVignettes("GenomicFeatures")

data(HepG2) # cell lines of liver origen
data("GM12878") # immortalized B cell
# Both being GRanges Objects


# Aim: to find overlaps in the two cell lines


HepG2_ind = queryHits(findOverlaps(HepG2, GM12878))
HepG2[HepG2_ind,] # showing the genes in HepG2 which overlaps that in GM12878


d_Hep_B = disjoin(c(HepG2,GM12878))
ind = which(d_Hep_B %over% HepG2 & d_Hep_B %over% GM12878)
d_Hep_B[ind] # showing the exact overlapping regions between the two cell lines

par(mfrow=c(3,1))
plotRanges(ranges(HepG2))
plotRanges(ranges(GM12878))
plotRanges(ranges(d_Hep_B[ind]))
# plot 11


# Finding Overlaps Assessment


# Where does the 17th HepG2 region start?

ranges(HepG2[17,])
#  [1]  46528596  46529164       569


# Closest regions in distinct GRanges

# Use distanceToNearest to find the closest region 
# in GM12878 to the 17th region in HepG2. 
# What is the start site of this region?

distanceToNearest(HepG2[17,], subject = GM12878)
#       queryHits subjectHits |  distance
#       <integer>   <integer> | <integer>
#   [1]         1         945 |      2284

ranges(GM12878[945,])
# [1]  46524762  46526311      1550


# Measuring distance between closest regions

# What is the distance between the 17th region of HepG2 
# and its closest region in GM12878?

mcols(distanceToNearest(HepG2[17,], subject = GM12878))$distance
# [1] 2284


# Summarizing proximities of nearest regions in a pair of GRanges

# For each region in HepG2 find the closest region in GM12878 
# and record the distance. 
# What proportion of these distances are smaller than 2000 base pairs? 
# Distance is a metadata column on the Hits object, so consider mcols().

d=mcols(distanceToNearest(HepG2, subject = GM12878))$distance
mean(d<2000)
# [1] 0.2673267


###


# Working with the SummarizedExperiment class


biocLite('airway')
library(airway)
# cells from human airway

data(airway)
airway
# class: RangedSummarizedExperiment 
# dim: 64102 8 
# metadata(1): ''
# assays(1): counts
# rownames(64102): ENSG00000000003 ENSG00000000005 ... LRG_98 LRG_99
# rowData names(0):
# colnames(8): SRR1039508 SRR1039509 ... SRR1039520 SRR1039521
# colData names(9): SampleName cell ... Sample BioSample


# looking at the assays

assay(airway)[1:4,1:3]
#                 SRR1039508 SRR1039509 SRR1039512
# ENSG00000000003        679        448        873
# ENSG00000000005          0          0          0
# ENSG00000000419        467        515        621
# ENSG00000000457        260        211        263


# looking at the rowData

rowData(airway)
# DataFrame with 64102 rows and 0 columns


# looking at the colData

names(colData(airway))
# [1] "SampleName" "cell"       "dex"        "albut"      "Run"        "avgLength" 
# [7] "Experiment" "Sample"     "BioSample"

# These are actually phenoData

colData(airway)[1:4,]
#            SampleName     cell      dex    albut        Run avgLength Experiment
#              <factor> <factor> <factor> <factor>   <factor> <integer>   <factor>
# SRR1039508 GSM1275862   N61311    untrt    untrt SRR1039508       126  SRX384345
# SRR1039509 GSM1275863   N61311      trt    untrt SRR1039509       126  SRX384346
# SRR1039512 GSM1275866  N052611    untrt    untrt SRR1039512       126  SRX384349
# SRR1039513 GSM1275867  N052611      trt    untrt SRR1039513        87  SRX384350
#                Sample    BioSample
#              <factor>     <factor>
#  SRR1039508 SRS508568 SAMN02422669
#  SRR1039509 SRS508567 SAMN02422675
#  SRR1039512 SRS508571 SAMN02422678
#  SRR1039513 SRS508572 SAMN02422670

# Note the paired design for the $dex
# i.e. untrt vs. trt


# Note that this is a RangedSummarizedExperiment

# We can get the ranges for the features
# using rowRanges()

gr_airway = rowRanges(airway)
class(gr_airway)
# [1] "GRangesList"
# attr(,"package")
# [1] "GenomicRanges"

gr_airway


# Considering the gene ORMDL3 with the ENSEMBL identifier ENSG00000172057

gr_ORMDL3 = rowRanges(airway)$ENSG00000172057
summary(gr_ORMDL3)
# [1] "GRanges object with 20 ranges and 2 metadata columns"

# Finding the overlapping regions using reduce()

summary(reduce(gr_ORMDL3))
# [1] "GRanges object with 8 ranges and 0 metadata columns"


###


# Using SummarizedExperiment with 450k Methylation arrays


biocLite('IlluminaHumanMethylation450kmanifest')
biocLite('IlluminaHumanMethylation450kanno.ilmn12.hg19')


biocLite('ArrayExpress')
library(ArrayExpress)

getAE('E-MTAB-5797')

dir(pattern = '5797') 
# [1] "E-MTAB-5797.idf.txt"   "E-MTAB-5797.raw.1.zip" "E-MTAB-5797.sdrf.txt" 
dir(pattern = 'idat')
# the idat are the min. processed data


# Reading the *.txt

library(data.table)
sd5797 = fread("E-MTAB-5797.sdrf.txt")
head(sd5797[,c(3,16,18)])
#     Characteristics[age] Label                   Performer
#  1:                   58   Cy3 IntegraGen SA, Evry, France
#  2:                   58   Cy5 IntegraGen SA, Evry, France
#  3:                   72   Cy3 IntegraGen SA, Evry, France
#  4:                   72   Cy5 IntegraGen SA, Evry, France
#  5:                   70   Cy3 IntegraGen SA, Evry, France
#  6:                   70   Cy5 IntegraGen SA, Evry, France


# And you can get the *.idat related by

sd5797$`Array Data File`
#  [1] "9406922003_R04C02_Grn.idat" "9406922003_R04C02_Red.idat"
#  [3] "9406922003_R02C01_Grn.idat" "9406922003_R02C01_Red.idat"
#  [5] "9406922003_R05C01_Grn.idat" "9406922003_R05C01_Red.idat"
#  [7] "9406922003_R01C01_Grn.idat" "9406922003_R01C01_Red.idat"
#  [9] "9406922003_R03C02_Grn.idat" "9406922003_R03C02_Red.idat"

# So that you would not mix with other files when you have getAE() for multiple times


# Using gsub() to get the prefix

pref = unique(gsub('_Grn.idat','',sd5797$`Array Data File`))
pref = unique(gsub('_Red.idat','',pref))
pref


# Or using substr()

pref_2 = unique(substr(sd5797$`Array Data File`,1,17))
pref_2 == pref # [1] TRUE TRUE TRUE TRUE TRUE


# Loading the data using the 'minfi' package and the read.metharray() command

library(BiocInstaller)
biocLite('minfi')
library(minfi)

raw = read.metharray(pref)
raw
# class: RGChannelSet 
# dim: 622399 5 
# metadata(0):
# assays(2): Green Red
# rownames(622399): 10600313 10600322 ... 74810490 74810492
# rowData names(0):
# colnames(5): 9406922003_R04C02 9406922003_R02C01 9406922003_R05C01
#   9406922003_R01C01 9406922003_R03C02
# colData names(0):
# Annotation
#   array: IlluminaHumanMethylation450k
#   annotation: ilmn12.hg19


# Generate SummarizedExperiment using preprocessQuantile()
# also from the minfi package
# convert the raw to normalized data


glioMeth = preprocessQuantile(raw)
glioMeth
# class: GenomicRatioSet 
# dim: 485512 5 
# metadata(0):
# assays(2): M CN
# rownames(485512): cg13869341 cg14008030 ... cg08265308 cg14273923
# rowData names(0):
# colnames(5): 9406922003_R04C02 9406922003_R02C01 9406922003_R05C01
#   9406922003_R01C01 9406922003_R03C02
# colData names(1): predictedSex
# Annotation
#   array: IlluminaHumanMethylation450k
#   annotation: ilmn12.hg19
# Preprocessing
#   Method: Raw (no normalization or bg correction)
#   minfi version: 1.24.0
#   Manifest version: 0.4.0

# Note that there are 2 assays
# One is M: methylation
# CN: copy number

assay(glioMeth[1:4,1:3], 'M')
#            9406922003_R04C02 9406922003_R02C01 9406922003_R05C01
# cg13869341         1.7581740         1.6586446         1.8828883
# cg14008030         0.7279845         0.4959878         0.3320419
# cg12045430        -1.7021606        -1.6417207        -1.5100296
# cg20826792        -1.0123667        -1.0926002        -0.8586541

assay(glioMeth[1:4,1:3], 'CN')
#            9406922003_R04C02 9406922003_R02C01 9406922003_R05C01
# cg13869341          14.40168          14.39969          14.34198
# cg14008030          14.61342          14.59784          14.76475
# cg12045430          14.72333          14.71023          14.65453
# cg20826792          14.68389          14.57409          14.64232

# 'CN' being some sort of relative abundance measues
# that will have to be reinterpreted to understand the actual relationship
# of deletion / duplication


###


# DataFrame and SummarizedExperiment assessment


# Tabulating sample characteristics

# The erma package includes detailed information on cell lines 
# analyzed in the epigenomics road map project. 
# You can query anatomic locations from which samples were derived as follows

biocLite('erma')
library(erma)
ee = makeErmaSet()
class(colData(ee))
length(names(colData(ee)))  
# 95  # lots of attributes!
table(ee$ANATOMY)
#    BLOOD    BRAIN      ESC      FAT     IPSC    LIVER     LUNG     SKIN VASCULAR 
#       15        7        2        1        1        1        2        1        1


# How many samples are derived from brain?
# 7


# DataFrame columns of arbitrary type

# Use the ErmaSet instance generated in the previous problem. Consider the code

mydf = colData(ee)[,1:10]
getClass("DataFrame")
mydf$demomat = matrix(0, nr=nrow(mydf), nc=5)
dim(mydf$demomat)
dim(mydf)
dim(data.frame(mydf))

# Why do the last two dim() results disagree?

# The data.frame class cannot treat the matrix 'demomat' 
# as an atomic item but breaks it by columns into a list; 
# each list entry counts as a column in the data.frame


###

# Memory-sparing approaches with HDF5 and indexed files

###


# External HDF5


library(BiocInstaller)
biocLite('HDF5Array')
library(HDF5Array)
biocLite('airway')
library(airway)


data("airway")
airway
Airway_assay = assay(airway)


# Writing HDF5Array

writeHDF5Array(Airway_assay,'aw.h5','airway')
file.exists('aw.h5') # [1] TRUE


# Reading HDF5Array

HDF5Array('aw.h5','airway')
# HDF5Matrix object of 64102 x 8 integers:
# ...


# Saving the whole SummarizedExperiment:

saveHDF5SummarizedExperiment(airway,'ext_airway',replace = T)
ext_airway = loadHDF5SummarizedExperiment('ext_airway')
# So we are having an external source of the SummarizedExperiment


# The delayed assay

assay(ext_airway)
# DelayedMatrix object of 64102 x 8 integers:
# ...


# That means it would not be calculated now, but it will when called


###


# GenomicFiles: application to many BAM files

# Important function: the MAP function


biocLite('GenomicFiles')
library(GenomicFiles)


biocLite('RNAseqData.HNRNPC.bam.chr14')
library(RNAseqData.HNRNPC.bam.chr14)


gf = GenomicFiles(files = RNAseqData.HNRNPC.bam.chr14_BAMFILES)
files(gf)  # getting all the files directory


# Binding the region of interest to the object

hn = GRanges("chr14", IRanges(21677296, 21737638), strand="-")
rowRanges(gf) = hn


rowData(gf) # DataFrame with 1 row and 0 columns
colData(gf) # DataFrame with 8 rows and 0 columns
# still nothing down there


library(GenomicAlignments)
MAP = function(r, f) {
  readGAlignmentPairs(f, param=ScanBamParam(which=r))
}
# r referring to the range of interest, 
# and f referring to the file being parsed for alignments overlapping the range.
ali = reduceByRange(gf, MAP=MAP)
ali

sapply(ali[[1]], length)
# ERR127306 ERR127307 ERR127308 ERR127309 ERR127302 ERR127303 ERR127304 ERR127305 
#      2711      3160      2946      2779        86        98       158       141

# So the first four samples are the wildtype samples
# And the latter four samples are the knockout samples


###


# Multiple BED files: a slice of the epigenomics roadmap


# The book chapter associated with this screencast discusses a 
# mild extension of GenomicFiles class to allow definition of some special methods.  
# The show method for instance includes some hints for how to work 
# with the packaged slice of the Epigenomics Road Map data.  
# For further details on some of the annotations and interpretations 
# of data visualized at the end, see section 2.2 of the erma vignette.


library(erma)

# The erma is developed as a demonstration of the utility of 
# packaging voluminous data on epigenomic assay outputs 
# on diverse cell types for 'out of memory' analysis.

# https://www.bioconductor.org/packages/release/bioc/vignettes/erma/inst/doc/erma.html


erset = makeErmaSet()
erset
# ErmaSet object with 0 ranges and 31 files:
# ...

colData(erset)
# (showing narrow slice of  31 x 95  DataFrame)   narrDF with 31 rows and 6 columns
# ...

table(erset$ANATOMY)
#    BLOOD    BRAIN      ESC      FAT     IPSC    LIVER     LUNG     SKIN VASCULAR 
#       15        7        2        1        1        1        2        1        1


# Using the stateProfile() function
# Suppose we are going to look at the range: 26:31
erset$Standardized.Epigenome.name[26:31]
# [1] "Brain Germinal Matrix"                "Brain Hippocampus Middle"            
# [3] "Brain Inferior Temporal Lobe"         "Brain_Dorsolateral_Prefrontal_Cortex"
# [5] "Fetal Lung"                           "Lung"     

stateProfile(erset[,26:31], shortCellType=FALSE)

# Plot 12
# The black line would be the starting point of transcription
# and to their left would be the promotor region, and 
# we can see different achitecture across different cell type

# Also note the [+] sign on the R to IL33
# that means a positive strand, and the transcription should be from left to right

stateProfile(erset, 'IL33', shortCellType = F)

# plot 13
# just typing in the gene 'IL33'


# Looking deeping into the 'state' in the stateProfile
# Suppose we are looking at the promotor: 14_EnhA2

data("states_25")
states_25[14,]
#    STATENO. MNEMONIC       DESCRIPTION COLOR.NAME COLOR.CODE     rgb
# 14       14    EnhA2 Active Enhancer 2     Orange 255,195,77 #FEC24D

# And we can have the description of it


###


# DNA Variants in 2000+ genomes: VcfStack class


library(BiocInstaller)
biocLite('ldblock')
library(ldblock)


# Loading the 1000 genomes project

sta = stack1kg()
sta
# ?failed here


# add information about geographic origin of samples

library(ph525x)
data(sampInfo_1kg)
rownames(sampInfo_1kg) = sampInfo_1kg[,1] # the HGxxxxx
cd = sampInfo_1kg[ colnames(sta), ]
colData(sta) = DataFrame(cd)


# Importing and using VCF data

library(erma)
orm = range(genemodel("ORMDL3")) # quick lookup

genome(orm) = "b37"  # must agree with VCF
seqlevelsStyle(orm) = "NCBI"  # use integer chromosome codes
ormRead = readVcfStack(sta, orm)
ormRead
# class: CollapsedVCF 
# dim: 179 2504
# ..

library(VariantTools)
vr = as(ormRead[,1:5], "VRanges")
vr[1:3,]
# VRanges object with 3 ranges and 28 metadata columns:
#               seqnames               ranges strand         ref
#                  <Rle>            <IRanges>  <Rle> <character>
#   rs145809834       17 [38077298, 38077298]      *           C
#   rs141712028       17 [38077366, 38077366]      *           C
#     rs3169572       17 [38077412, 38077412]      *           G
# ...

table(vr[which(sampleNames(vr)=="HG00096"),]$GT)
# 0|0 1|1 
# 175   4


###


# External data resources assessment


# Using HDF5 with a SummarizedExperiment

# Let's use the airway package again to obtain a SummarizedExperiment, 
# saving it as HDF5 in a temporary location.

library(airway)
td = tempfile()
data("airway")
saveHDF5SummarizedExperiment(airway, td)

# After this save operation completes, what is

length(dir(td)) # [1] 2
dir(td) # [1] "assays.h5" "se.rds" 


# Continuation: components of the saved SummarizedExperiment

X  = readRDS(dir(td, full=TRUE)[2]) 

# What is the class of X?

class(X) # RangedSummarizedExperiment
X
# class: RangedSummarizedExperiment 
# dim: 64102 8 
# ..
rowData(X)
# DataFrame with 64102 rows and 0 columns


# Hazards of bypassing the standard approach to reloading

# At this stage, X looks like a familiar object.
# What is the class of assay(X[1,1])?

# try-error


# Counting RNA-seq reads over an exon

# In our illustration of RNA-seq read tabulation in the HNRNPC knockdown experiment, 
# we used the whole coding region of HNRNPC for counting. 
# Instead, we can use the gene model to find exons and focus our attention on these. 
# The following code can be used for this:

library(erma)
hn = genemodel("HNRNPC")
hn
# GRanges object with 19 ranges and 2 metadata columns:
#        seqnames               ranges strand |   exon_id      SYMBOL
# ...

e1 = hn[1]  # first exon in the model

library(GenomicFiles)
library(RNAseqData.HNRNPC.bam.chr14)
gf = GenomicFiles(files=RNAseqData.HNRNPC.bam.chr14_BAMFILES)
gf
rowRanges(gf) = e1 # And this would be the range of interest

library(GenomicAlignments)
MAP = function(r, f) {
  readGAlignmentPairs(f, param=ScanBamParam(which=r))
}
ali = reduceByRange(gf, MAP=MAP)
elementNROWS(ali[[1]]) # use [[1]] as there is only one request
# ERR127306 ERR127307 ERR127308 ERR127309 ERR127302 ERR127303 ERR127304 ERR127305 
#       497       617       482       569         6         3         9        23

# What is the count of reads found to align over the first exon in the first sample (ERR127306)?


# Querying BED files to compare epigenomic states of cell types

# Use the code below to generate a visualization of regulatory states 
# of chromatin in the vicinity of the transcriptional start site (TSS) 
# of gene CD28.

library(erma)
ermaset = makeErmaSet()
stateProfile(ermaset[,c(4,6,30,31)], "CD28", short=FALSE)
# c(4,6,30,31) specifying the sample sites (fetal lung, lung, primary B and T cells)

# plot 14
# Primary T and B cells obtained from peripheral blood are 
# discordant with respect to promoter or enhancer state upstream of CD28


###

# Multi-omics solutions; role of cloud resources  

###


# Collecting multiple molecular assay outputs on a set of samples


biocLite('MultiAssayExperiment')
biocLite('RaggedExperiment')

library(MultiAssayExperiment)
library(RaggedExperiment)


# The download.file() has some sort of error
# we therefore need to manually download the file
# "http://s3.amazonaws.com/multiassayexperiments/gbmMAEO.rds"

gbm = readRDS("gbmMAEO.rds")
library(SummarizedExperiment)
gbm = updateObject(gbm)
gbm


biocLite('UpSetR')
library(UpSetR)
upsetSamples(gbm)
# plot 15
# Set size and Intersection size


# As the .rds file has some error, the following codes were just copies from the course

mut = experiments(gbm)[["Mutations"]]
mut
## class: RaggedExperiment 
## dim: 22073 290 
## assays(73): Entrez_Gene_Id Center ... OREGANNO_ID OREGANNO_Values
## rownames(22073): ATAD3B TPM3 ... F9 FATE1
## colnames(290): TCGA-02-0003-01A-01D-1490-08
##   TCGA-02-0033-01A-01D-1490-08 ... TCGA-81-5911-01A-12D-1845-08
##   TCGA-87-5896-01A-01D-1696-08
## colData names(0):

head(assayNames(mut))
## [1] "Entrez_Gene_Id"         "Center"                
## [3] "NCBI_Build"             "Variant_Classification"
## [5] "Variant_Type"           "Reference_Allele"

rowRanges(mut)
## GRanges object with 22073 ranges and 0 metadata columns:
##          seqnames                 ranges strand
##             <Rle>              <IRanges>  <Rle>
##   ATAD3B        1 [  1430871,   1430871]      +
##     TPM3        1 [154148652, 154148652]      +
##    NR1I3        1 [161206281, 161206281]      +
##      AGT        1 [230846235, 230846235]      +
##    TACC2       10 [123810032, 123810032]      +
##      ...      ...                    ...    ...
##    MAGT1        X [ 77112989,  77112989]      +
##   RHOXF1        X [119243159, 119243159]      +
##    RAP2C        X [131348336, 131348336]      +
##       F9        X [138643011, 138643011]      +
##    FATE1        X [150891145, 150891145]      +
##   -------
##   seqinfo: 24 sequences from hg19 genome; no seqlengths

sort(table(names(rowRanges(mut))),decreasing=TRUE)[1:5]
## 
## Unknown     TTN    EGFR    TP53    PTEN 
##     126     121     102     101      93

table(as.character(assay(mut, "Variant_Classification")))
## 
##        Frame_Shift_Del        Frame_Shift_Ins           In_Frame_Del 
##                    566                    217                    214 
##           In_Frame_Ins      Missense_Mutation      Nonsense_Mutation 
##                     28                  14213                    851 
##       Nonstop_Mutation                 Silent            Splice_Site 
##                     17                   5514                    382 
## Translation_Start_Site 
##                     71

# Which genes have had deletions causing frame shifts?
# Given that assay(mut, "Variant_Classification") is a matrix, we can use apply over rows

rfs = rowRanges(mut)[ which( apply(
  assay(mut, "Variant_Classification"), 1, 
  function(x) any(x=="Frame_Shift_Del"))
)]
rfs
## GRanges object with 566 ranges and 0 metadata columns:
##            seqnames                 ranges strand
##               <Rle>              <IRanges>  <Rle>
##    ZNF280D       15   [56993158, 56993158]      +
##     CREBBP       16   [ 3843446,  3843446]      +
##     SPACA3       17   [31322643, 31322643]      +
##     GGNBP2       17   [34943625, 34943625]      +
##      MALT1       18   [56400716, 56400716]      +
##        ...      ...                    ...    ...
##    TBC1D8B        X [106066520, 106066521]      +
##     SPTBN5       15 [ 42164092,  42164092]      +
##   C16orf82       16 [ 27078770,  27078770]      +
##   C12orf42       12 [103695960, 103695960]      +
##      DUSP8       11 [  1577819,   1577820]      +
##   -------
##   seqinfo: 24 sequences from hg19 genome; no seqlengths


###


# MultiAssayExperiment
# exemplified by another dataset
## also Multiomic TCGA data assessment


laml = readRDS("tcgaLAML.rds")
laml
# A MultiAssayExperiment object of 5 listed
#  experiments with user-defined names and respective classes. 
#  Containing an ExperimentList class object of length 5: 
# [1] RNASeqGene: ExpressionSet with 19990 rows and 179 columns 
# [2] RNASeq2GeneNorm: ExpressionSet with 20501 rows and 173 columns 
# [3] CNASNP: RaggedExperiment with 874897 rows and 392 columns 
# [4] CNVSNP: RaggedExperiment with 28324 rows and 380 columns 
# [5] Methylation: SummarizedExperiment with 27578 rows and 194 columns 
# Features: 

biocLite('UpSetR')
library(UpSetR)
upsetSamples(laml)
# plot 16

# What is the number of samples providing CNASNP, CNVSNP, and methylation, but no RNA-seq data?
# 17


# using the experiments()
# What is
length(experiments(laml))
# [1] 5


# supposedly we would like to look into the Methylation
experiments(laml)
met = experiments(laml)[['Methylation']]
met
# class: SummarizedExperiment 
# dim: 27578 194 
# metadata(0):
# assays(1): ''
# rownames(27578): cg00000292 cg00002426 ... cg27662877 cg27665659
# rowData names(3): Gene_Symbol Chromosome Genomic_Coordinate
# colnames(194): TCGA-AB-2802-03A-01D-0741-05 TCGA-AB-2803-03A-01D-0741-05 ... TCGA-AB-3011-03A-01D-0742-05
#   TCGA-AB-3012-03A-01D-0741-05
# colData names(0):


CNASNP = experiments(laml)[['CNASNP']]
CNASNP
assayNames(CNASNP)
# [1] "Num_Probes"   "Segment_Mean"
rR_CNASNP = rowRanges(CNASNP)
rR_CNASNP
# GRanges object with 874897 ranges and 0 metadata columns:
#                         seqnames              ranges strand
#                            <Rle>           <IRanges>  <Rle>
#         1:61735-804456        1        61735-804456      *
#       1:809773-1478153        1      809773-1478153      *
# ...


# Ranges in the assay
sort(table(names(rowRanges(CNASNP))),decreasing=TRUE)[1:5]
#  13:72477555-72480556 5:103860373-103860604   3:68746259-68747413  15:20581451-20581564  21:44970873-44972812 
#                   172                   153                   137                   125                   114


# Tabulating the assays
sort(table(as.character(assay(CNASNP, 'Num_Probes'))),decreasing = T)[1:10]
#     8     7     6     9    10     5    11    12     4    13 
# 33657 33616 32522 31959 30377 29974 28547 26408 26090 24043


# Quick check for aberrant samples

# We will use principal components analysis to assess the normalized RNA-seq data for sample aberrations.

# The following code obtains a log-transformed version of the normalized gene-level RNASeq data, 
# and uses all measured genes to re-express sample-to-sample variation through principal components. 
# We use a pairs plot of five PCs to look for clumping of samples in 2-dimensional projections.

lnorna = log(exprs(experiments(laml)[["RNASeq2GeneNorm"]])+1)
pca = prcomp(t(lnorna))
pairs(pca$x[,1:5])

# How many samples are separated from the bulk of the data by having high values of PC5?
# 4


# looking deeper
experiments(laml)[["RNASeq2GeneNorm"]]
dim(pca$x) # [1] 173 173
dim(pca$rotation) # [1] 20501   173
# So the clusterings were between samples


# Isolating another unusual group

# We can also see some clumping of samples when PC3 is plotted against PC4. 
# One way to identify the samples is to use residuals from robust regression. 
# We'll use the lqs procedure from the MASS package 
# that uses approximate minimization of the median of squared residuals 
# to obtain robust estimates of regression coefficients. 
# The idea is that the procedure finds the line that best fits the "majority" of points. 
# Since the clump we are interested in is clearly a minority, 
# these points will not play a substantial role in determining slope of the lqs fit, 
# and will have relatively large residuals.

library(MASS)
set.seed(1234)
zz2 = lqs(PC3~PC4, data=data.frame(pca$x))

# Now we set up a multipanel display, and plot the PCs, the fitted line, the residuals, 
# and the isolated aberrant points.

par(mfrow=c(2,2), mar=c(4,4,1,1))
plot(PC3~PC4, data=data.frame(pca$x), main="PC3 vs PC4")
abline(coef(zz2), col="green", lwd=2)
legend(-40,64, lty=1, col="green", legend="lqs fit", bty="n")
qqnorm(resid(zz2), main="lqs residuals")
abline(h= -50)
plot(PC3~PC4, data=data.frame(pca$x), main="isolate by resid < -50 and PC4>0")
points(PC3~PC4, data=data.frame(pca$x),
       subset=which(resid(zz2) < -50 & pca$x[,"PC4"] > 0), pch=19, col="red")
# plot 17

# How many samples are identified in this procedure?
# 16


###

# Section 3: Genomic annotation with Bioconductor

# Prologue: Detailed applications of GRanges

###


# Setup and GRanges check


# This week we are using two sets of genomic regions as examples. 
# These are the reported Estrogen Related Receptor binding sites obtained 
# for a ENCODE ChIP-seq experiment on two cell lines: HepG2 and GM12787. 
# We have put these regions into an R package for your convenience. 
# If you have not done so already, please download and install the package:

library(BiocInstaller)
biocLite("genomicsclass/ERBS")
library(ERBS)


# Setup/version check
# checking up the latest version 

library(BiocInstaller)
biocVersion()
# [1] '3.6'


# Classes for genomic ranges

# Load the ERBS library and then the HepG2 object.
# What is the class of HepG2 (hint: use the class function)?

data(HepG2) 
class(HepG2)
# [1] "GRanges"  # class of HepG2
# attr(,"package")
# [1] "GenomicRanges"  # class of package HepG2 is defined in


# Counting regions

# Explore the HepG2 object.
# How many regions are represented?

HepG2
# GRanges object with 303 ranges and 7 metadata columns:
# ...


###


# Introduction to Genomic Ranges Assessment


# Statistics on peak scores

# What is the median of the signalValue column for the HepG2 data?

median(mcols(HepG2)$signalValue)
# [1] 7.024


# Locating the largest peak

# In what chromosome is the region with the highest signalValue (copy and paste your answer)?

ind = which(rank(-mcols(HepG2)$signalValue)==1)  # note the '-' before the mcols()
seqnames(HepG2)[ind]
# Values : chrX
HepG2$signalValue[ind]
# [1] 91.779


# Tabulating by chromosome

# How many regions are from chromosome 16?

ind = which(seqnames(HepG2)=='chr16')
length(ind)  # [1] 31


# Statistics on peak widths

# Make a histogram of the widths of the regions from all chromosomes (not just chr16). 
# Note it has a heavy right tail.

HepG2_width = width(HepG2)
hist(HepG2_width)
# plot 18

# What is the median width?

median(HepG2_width) # [1] 560


# Other techniques that are taught by the intro video


# 'Rle' to string

seqnames(HepG2[1,])
# factor-Rle of length 1 with 1 run
# Lengths:    1
# Values : chr2
# Levels(93): chr1 chr2 chr3 chr4 chr5 chr6 ... chrUn_gl000245 chrUn_gl000246 chrUn_gl000247 chrUn_gl000248 chrUn_gl000249

# Note that it would become something like this

as.character(seqnames(HepG2[1,]))
# [1] "chr2"
# This would be better

# Also changing the Rle to chr can give us a nicer table
# or else some figures e.g. chrUn_gl000212 would be present
chr = as.character(seqnames(HepG2))
table(chr)
#  chr1 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19  chr2 chr20 chr21 chr22  chr3  chr4  chr5  chr6  chr7 
#    18     9     9    13     2     8     5    31    21     6    16    38    27     3     5    15     4     8    24    14 
#  chr8  chr9  chrX 
#    11    12     4


# [order()] in GRanges objects

HepG2[order(HepG2)]
# And the it would be ordered according to the seqnames and ranges

# And the much better arrangements for the seqnames

x = HepG2[order(HepG2)]
seqnames(x)
# factor-Rle of length 303 with 23 runs
#   Lengths:    18    38    15     4     8    24    14    11    12 ...    31    21     6    16    27     3     5     4
#   Values :  chr1  chr2  chr3  chr4  chr5  chr6  chr7  chr8  chr9 ... chr16 chr17 chr18 chr19 chr20 chr21 chr22  chrX
# Levels(93): chr1 chr2 chr3 chr4 chr5 chr6 ... chrUn_gl000245 chrUn_gl000246 chrUn_gl000247 chrUn_gl000248 chrUn_gl000249


###


# Genes as GRanges Assessment


biocLite('Homo.sapiens')
library(Homo.sapiens)
ghs = genes(Homo.sapiens)  # Note the use of the genes() command


ghs
#         seqnames              ranges strand |       GENEID
#            <Rle>           <IRanges>  <Rle> | <FactorList>
#       1    chr19   58858172-58874214      - |            1
#      10     chr8   18248755-18258723      + |           10
# ...


# What genome build was used here?
Homo.sapiens
# ...
# # Based on genome:  hg19 
# ...

# Or you can simply use
genome(ghs)


# How many genes are represented in ghs?

ghs
# GRanges object with 23056 ranges and 1 metadata column:
# ...

# or
length(ghs) # [1] 23056


# Sorting genes into chromosomes

# What is the chromosome with the most genes?

chr = as.character(seqnames(ghs))
table(chr)
#                  chr1                 chr10                 chr11                 chr12                 chr13 
#                  2326                   903                  1439                  1173                   449
# ....
which.max(table(chr))
# chr1 
# 1


# The distribution of "gene lengths"

# Make a histogram of the widths of genes (use the width() on the GRanges object). 
# This width gives the number of basepairs from the start of the gene to the end, so including exons and introns.
# Which best describes the width of genes?

hist(width(ghs))
# plot 19
# A distribution skewed to the right


# Statistics on gene lengths

# What is the median gene width?

median(width(ghs))
# [1] 20115.5


# Other techniques in the intro video


# precede()

library(ERBS)
data(HepG2) 
HepG2_GR = GRanges(HepG2)

precede(ghs)[1]  # [1] 19603
ranges(ghs[precede(ghs)[1],])  # 79673  58637619  58666477     28859
ranges(ghs[1])  # 1  58858172  58874214     16043
# This can show the ranges which precedes the respected range


ind = precede(ghs,subject = HepG2_GR)
ghs[1]
# GRanges object with 1 range and 1 metadata column:
#   seqnames            ranges strand |       GENEID
#      <Rle>         <IRanges>  <Rle> | <FactorList>
# 1    chr19 58858172-58874214      - |            1
HepG2_GR[ind[1],]
# GRanges object with 1 range and 7 metadata columns:
#     seqnames            ranges strand |      name     score       col signalValue    pValue       qValue      peak
#        <Rle>         <IRanges>  <Rle> | <numeric> <integer> <logical>   <numeric> <numeric>    <numeric> <integer>
# [1]    chr19 55972235-55973446      * |      <NA>         0      <NA>       5.428     21.58 2.613574e-19       411

# ghs[1] precedes subject[ind[1]]


###


#  Finding and getting annotation for closest gene Assessment


library(ERBS)
data(HepG2)
data(GM12878)
res = findOverlaps(HepG2,GM12878)
erbs = HepG2[queryHits(res)]
erbs = granges(erbs)

# The following command is similar:
  
erbs2 = intersect(HepG2,GM12878)


# Comparing consensus methods

# Which of the following is true of the consensus entities erbs and erbs2?

x = erbs[order(erbs)]
y = erbs2[order(erbs2)]

x==y
mean(x==y)  # [1] 0.9333333

par(mfrow=c(2,1))
library(ph525x)
plotRanges(x)
plotRanges(y)
# plot 20

# Over 90% of these regions in these two objects are the same with the different regions being smaller in erbs2.

# Or we can use this code

mean(start(x)==start(y) & end(x)==end(y))
# [1] 0.9333333
mean(width(x)>=width(y))
# [1] 1


# A one-liner for transcription start sites
## Introducing the function resize()

# Using the ghs regions:
  
library(Homo.sapiens)
ghs = genes(Homo.sapiens)

# and what you learned in the video, convert the ghs object to one that represents just the tss
# What is the TSS (Transcription Start Site) of the gene with ID: 100113402?

# Hint: look at the ghs in the console. 
# Note that the names of the ranges are the same as the GENEID column. 
# So you can index the ranges directly with "100113402"

geneid = as.data.frame(ghs$GENEID)
ind = which(geneid$value=='100113402')
ghs[ind,]
# GRanges object with 1 range and 1 metadata column:
#           seqnames            ranges strand |       GENEID
#              <Rle>         <IRanges>  <Rle> | <FactorList>
# 100113402    chr16 70563402-70563502      + |    100113402

# simply:
ranges(resize(ghs[ind,],1))


# The gene with TSS nearest a binding site
## Using the distanceToNearest()

# Now using the erbs regions defined in a previous question:
  
library(ERBS)
data(HepG2)
data(GM12878)
res = findOverlaps(HepG2,GM12878)
erbs = HepG2[queryHits(res)]
erbs = granges(erbs)

# What is the GENEID of the gene with TSS closest to the 4th region of erbs?

tss = resize(ghs,1)
d = distanceToNearest(erbs,subject = tss)
ind = subjectHits(d[4])
ghs$GENEID[[ind,]]
# [1] 2101

# Extracting the values from d

queryHits(d)
subjectHits(d)
values(d)$distance


# In the question above, you identified a gene.
# Use the select function to determine which is the SYMBOL of this gene
## Using keytypes() and columns() to access other aspects of the database

keytypes(Homo.sapiens)
# ...  # so you can see what would be the possible keys
columns(Homo.sapiens)
# ... # And here you can see the 'columns' required in the select function
key = as.character(values(ghs)$GENEID)
keytype = 'GENEID'
col = 'SYMBOL'

sel_ghs = select(Homo.sapiens, keys = key, columns = col, keytype = keytype)
sel_ghs[which(sel_ghs$GENEID=='2101'),]
#      GENEID SYMBOL
# 6316   2101  ESRRA


# Other techniques from the intro video:


# Finding the genes with d<1000

ind = subjectHits(d)[which(values(d)$distance<=1000)]  # very important
# Need to get the index in the ghs set but not the erbs set
sel_ghs = select(Homo.sapiens,
                 keys = as.character(values(ghs)$GENEID),
                 keytype = 'GENEID',
                 columns = c('SYMBOL','GENENAME'))
sel_ghs[ind,]


###


# Getting Sequence Assessment


library(ERBS)
library(GenomicRanges)
data(HepG2)
data(GM12878)
res = findOverlaps(HepG2,GM12878)
erbs = HepG2[queryHits(res)]
erbs = granges(erbs)


# GC content of binding regions

# Now load the human genome data

biocLite('BSgenome.Hsapiens.UCSC.hg19')
library(BSgenome.Hsapiens.UCSC.hg19)
Hsapiens
class(Hsapiens)
# [1] "BSgenome"
# attr(,"package")
# [1] "BSgenome"

# Now use the getSeq function to extract the sequence of each region in erbs. 
# Then compute the GC-content 
# (the number of C's + the number of G's divided by the length of sequence) of each.

# Q:What is the median GC-content

?getSeq
# x	  A BSgenome object or any other supported object. Do showMethods("getSeq") to get the list of all supported types for x.
# ...	Any additional arguments needed by the specialized methods.

showMethods(getSeq)
# x="BSgenome"
# x="FaFile"
# x="FaFileList"
# x="TwoBitFile"
# x="XStringSet"

erbs_seq = getSeq(Hsapiens, erbs)
erbs_seq
# A DNAStringSet instance of length 75
# ...

# Calculating the GC-content
C = vcountPattern('C',erbs_seq)  # no need to use apply() or sapply()  # also note the 'v' before countPattern as there are multi strings
G = vcountPattern('G',erbs_seq)
CG_content = (C+G) / width(erbs_seq)
median(CG_content) # [1] 0.652568

# alternative method
alpha = alphabetFrequency(erbs_seq)
CG = (alpha[,2]+alpha[,3]) / width(erbs_seq)
median(CG) # [1] 0.652568
  

# GC content of a shifted set of regions

# Now create a control set of regions by shifting erbs by 10000.
# What is the median GC-content of these control regions:

erbs_c = shift(erbs,10000)
erbs_c_seq = getSeq(Hsapiens,erbs_c)
alpha = alphabetFrequency(erbs_c_seq)
CG = (alpha[,2]+alpha[,3]) / width(erbs_c_seq)
median(CG) # [1] 0.4860174


# Other techniques in the intro vid


# Matching a specific pattern

pat = 'TCAAGGTCA'
sum(vcountPattern(pat,erbs_seq)) # [1] 8

# Matching a reverse complement

sum(vcountPattern(pat,reverseComplement(erbs_seq))) # [1] 21


###


# Annotation of genes and transcripts 


###


# Reference genomes: assessment


# Reference genome discovery

# How many Bioconductor packages provide reference genomic sequence for zebrafish? 
# Exclude the packages with suffix .masked, that we will discuss later.:

library(Biostrings)
library(BSgenome)
ag = available.genomes()  # From the package 'BSgenome'
length(ag)

a = grep('Drerio',ag)
b = grep('.masked',ag,invert = T)
ag[a[a %in% b]]
# [1] "BSgenome.Drerio.UCSC.danRer10" "BSgenome.Drerio.UCSC.danRer5"  "BSgenome.Drerio.UCSC.danRer6" 
# [4] "BSgenome.Drerio.UCSC.danRer7"

# A better method:
grep("mask", grep("Drerio", available.genomes(), value=TRUE), invert=TRUE, value=TRUE) # exclude masked


# Masking structures for genome gaps and repetitions

# We have noted that the reference genome builds for complex organisms are works in progress. 
# Genomic sequence "mask" structures have been defined to isolate ambiguous, unmappable, 
# and low-complexity segments of genomes so that sequence analysis research can be targeted 
# to reflect current knowledge of sequence regions that are more likely to be functionally informative.

# Obtain BSgenome.Hsapiens.UCSC.hg19.masked (it is only a 20MB transfer.)

library(BiocInstaller)
biocLite('BSgenome.Hsapiens.UCSC.hg19.masked')

library(BSgenome.Hsapiens.UCSC.hg19.masked)
c17m = BSgenome.Hsapiens.UCSC.hg19.masked$chr17

# What is the class of c17m?
c17m
class(c17m)
# [1] "MaskedDNAString"
# attr(,"package")
# [1] "Biostrings"


# Quantifying assembly gaps

c17m
#   81195210-letter "MaskedDNAString" instance (# for masking)
# seq: AAGCTTCTCACCCTGTTCCTGCATAGATAATTGCATGACAATTGCCTTGTCCCTGCT...GGTTAGGGTGTGGGTGTGGGTGTGGGTGTGGGTGTGGTGTGTGGGTGTGGGTGTGGT
# masks:
#   maskedwidth maskedratio active names                               desc
# 1     3400000  0.04187439   TRUE AGAPS                      assembly gaps
# ..

# In build hg19, what percentage of the length of chromosome 22 is occupied by "assembly gaps"? 
# Reply with an integer between 0 and 100.

c22m = BSgenome.Hsapiens.UCSC.hg19.masked[['chr22']]
c22m
sum(width(masks(c22m)$AGAPS))/length(c22m)
# [1] 0.3198546

masks(c22m)$AGAPS
# NormalIRanges object with 5 ranges and 0 metadata columns:
#         start       end     width
#     <integer> <integer> <integer>
# [1]         1  16050000  16050000
# [2]  16697851  16847850    150000
# [3]  20509432  20609431    100000
# [4]  50364778  50414777     50000
# [5]  51244567  51304566     60000


# Other techniques in the intro video


# The packages of reference genome for Human

grep('Hsap',ag,value=T)
#  [1] "BSgenome.Hsapiens.1000genomes.hs37d5" "BSgenome.Hsapiens.NCBI.GRCh38"       
#  [3] "BSgenome.Hsapiens.UCSC.hg17"          "BSgenome.Hsapiens.UCSC.hg17.masked"  
#  [5] "BSgenome.Hsapiens.UCSC.hg18"          "BSgenome.Hsapiens.UCSC.hg18.masked"  
#  [7] "BSgenome.Hsapiens.UCSC.hg19"          "BSgenome.Hsapiens.UCSC.hg19.masked"  
#  [9] "BSgenome.Hsapiens.UCSC.hg38"          "BSgenome.Hsapiens.UCSC.hg38.masked"


# 'methods' for the class

library(BSgenome.Hsapiens.UCSC.hg19)
class(Hsapiens) # [1] "BSgenome"
methods(class='BSgenome')
#  [1] $               [[              as.list         bsgenomeName    coerce         
#  [6] coerce<-        commonName      countPWM        export          getSeq         
# [11] injectSNPs      length          masknames       matchPWM        mseqnames      
# [16] names           organism        provider        providerVersion releaseDate    
# [21] releaseName     seqinfo         seqinfo<-       seqnames        seqnames<-     
# [26] show            snpcount        snplocs         SNPlocs_pkgname sourceUrl      
# [31] toString        vcountPattern   vcountPDict     Views           vmatchPattern  
# [36] vmatchPDict    


# The substr()

Hs_chrx = Hsapiens[['chrX']]
substr(Hs_chrx, 5e6, 5.1e6)
#   100001-letter "DNAString" instance
# seq: GCCTCAATGTCAGAATTATGCTGTTGCCCAAAATTGCTT...TCTACTAAAAATACAAAAATTAGCTGGGCATGGTGGTG


# Simple operations on sequences

nchar(Hs_chrx)
# [1] 155270560

alphabetFrequency(Hs_chrx)
#        A        C        G        T        M        R        W        S        Y 
# 45648952 29813353 29865831 45772424        0        0        0        0        0 
#        K        V        H        D        B        N        -        +        . 
#        0        0        0        0        0  4170000        0        0        0


# Summing up all the sequences

sum(as.numeric(unlist(lapply(1:24,function(x) {
  nchar(Hsapiens[[x]])
}))))
# [1] 3095677412
# And it took really long


# Using the parallel library

library(parallel)
detectCores() # [1] 4
options(mc.cores=4)

sum(as.numeric(unlist(mclapply(1:24,function(x) {
  nchar(Hsapiens[[x]])
}))))  # Using mclapply() instead of lapply()


###


# Gene and transcript model assessment


# Parsing a gene model display

library(ph525x)

# Use the modPlot function to visualize a model for ESR1,
# the estrogen receptor 1 gene.

modPlot("ESR1", useGeneSym=FALSE, collapse=FALSE) 
# plot 21

# Each of the linear structures (line segments with right-pointing arrow glyphs starting and ending with, 
# and sometimes interrupted by, little yellowish polygons) is a
# Transcript

# The linear structures are transcripts, 
# composed of untranslated regions at the ends, 
# exons (small yellow rectangles), and introns (arrows between exons). 
# The collection of transcripts makes up a "gene model".


# Enumerating transcipts

# Finding the number of transcripts in ESR1

biocLite('TxDb.Hsapiens.UCSC.hg19.knownGene')
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb = TxDb.Hsapiens.UCSC.hg19.knownGene

ghs = genes(Homo.sapiens)
key = as.character(values(ghs)$GENEID)
keytype = 'GENEID'
col = 'SYMBOL'

sel_ghs = select(Homo.sapiens, keys = key, columns = col, keytype = keytype)
sel_ghs[which(sel_ghs$SYMBOL=='ESR1'),]
#      GENEID SYMBOL
# 6312   2099   ESR1

length(transcripts(txdb, filter=list(gene_id="2099")))
# [1] 27


# Other techniques in the book


# Filtering from the .KnownGene using the exon() command

exons(txdb, columns=c("EXONID", "TXNAME", "GENEID"),
      filter=list(gene_id=c(100, 101)))
# GRanges object with 39 ranges and 3 metadata columns:
#      seqnames                 ranges strand |    EXONID
#         <Rle>              <IRanges>  <Rle> | <integer>
#  [1]    chr10 [135075920, 135076737]      - |    144421
# ...
#                                  TXNAME          GENEID
#                         <CharacterList> <CharacterList>
#    [1] uc009ybi.3,uc010qva.2,uc021qbe.1             101


# And we can make see all the columns by columns()

columns(txdb)
# [1] "CDSCHROM"   "CDSEND"     "CDSID"      "CDSNAME"    "CDSSTART"   "CDSSTRAND" 
# [7] "EXONCHROM"  "EXONEND"    "EXONID"     "EXONNAME"   "EXONRANK"   "EXONSTART" 
# [13] "EXONSTRAND" "GENEID"     "TXCHROM"    "TXEND"      "TXID"       "TXNAME"    
# [19] "TXSTART"    "TXSTRAND"   "TXTYPE" 


# The ENSEMBL annotation


library(ensembldb)
biocLite('EnsDb.Hsapiens.v75')
library(EnsDb.Hsapiens.v75)

edb = EnsDb.Hsapiens.v75
columns(edb)
#  [1] "ENTREZID"            "EXONID"              "EXONIDX"            
#  [4] "EXONSEQEND"          "EXONSEQSTART"        "GENEBIOTYPE"        
#  [7] "GENEID"              "GENENAME"            "GENESEQEND"         
# [10] "GENESEQSTART"        "INTERPROACCESSION"   "ISCIRCULAR"         
# [13] "PROTDOMEND"          "PROTDOMSTART"        "PROTEINDOMAINID"    
# [16] "PROTEINDOMAINSOURCE" "PROTEINID"           "PROTEINSEQUENCE"    
# [19] "SEQCOORDSYSTEM"      "SEQLENGTH"           "SEQNAME"            
# [22] "SEQSTRAND"           "SYMBOL"              "TXBIOTYPE"          
# [25] "TXCDSSEQEND"         "TXCDSSEQSTART"       "TXID"               
# [28] "TXNAME"              "TXSEQEND"            "TXSEQSTART"         
# [31] "UNIPROTDB"           "UNIPROTID"           "UNIPROTMAPPINGTYPE"
# So there are more columns indeed, including the protein, etc..


# filtering in the ensembldb

txs <- transcripts(edb, filter = GenenameFilter("ZBTB16"),
                   columns = c("protein_id", "uniprot_id", "tx_biotype"))
# for the genename 'ZBTB16'
txs
# GRanges object with 20 ranges and 5 metadata columns:
#                  seqnames                 ranges strand |      protein_id
#                     <Rle>              <IRanges>  <Rle> |     <character>
#  ENST00000335953       11 [113930315, 114121398]      + | ENSP00000338157
# ...
#                    uniprot_id              tx_biotype           tx_id   gene_name
#                   <character>             <character>     <character> <character>
#  ENST00000335953  ZBT16_HUMAN          protein_coding ENST00000335953      ZBTB16


# The completely unfiltered version:

txs = transcripts(edb)


# Accessing the columns
# e.g. "TXBIOTYPE"

txs$tx_biotype
#    [1] "processed_transcript"               "transcribed_unprocessed_pseudogene"
#    [3] "transcribed_unprocessed_pseudogene" "transcribed_unprocessed_pseudogene"
#    [5] "unprocessed_pseudogene"             "unprocessed_pseudogene"   
# ...

table(txs$tx_biotype)
#           3prime_overlapping_ncrna                          antisense 
#                                 29                              10058 
#                          IG_C_gene                    IG_C_pseudogene 
#                                 31                                 13 
#                          IG_D_gene                          IG_J_gene 
#                                 64                                 24 
#                    IG_J_pseudogene                          IG_V_gene 
#                                  6                                185 
#                    IG_V_pseudogene                            lincRNA 
#                                264                              12101
# ...


###


# Import/Export assessment


# miRNA target sites: pre-GRanges

# Bioconductor's rtracklayer package supports import and export of files in common genomic data formats. 
# The package includes a demonstration dataset of microRNA target sites.

library(rtracklayer)
data(targets)

# What is the class of targets?

targets
class(targets)
# [1] "data.frame"


# Checking essential metadata

# To what reference build do the chromosome, start, and end values in targets refer?

# it cannot be determined without forensic work
# Unfortunately it is not typical to add metadata to data.frame instances, 
# as one has to fill a column with relevant information.


# GRanges to bed

# We can create a GRanges instance from the targets data frame as follows

library(GenomicRanges)
head(targets)
#                  name          target chrom     start       end strand
# 555774     hsa-miR-16 ENST00000000412 chr12   8985197   8985217      -
# 415091 hsa-miR-509-3p ENST00000003084  chr7 117095440 117095461      +
# ...
mtar = with(targets,
            GRanges(chrom, IRanges(start,end), strand=strand,
                    targets=target, mirname=name))
mtar
# GRanges object with 2981 ranges and 2 metadata columns:
#       seqnames              ranges strand |         targets        mirname
#          <Rle>           <IRanges>  <Rle> |        <factor>       <factor>
#  [1]    chr12     8985197-8985217      - | ENST00000000412     hsa-miR-16
# ...

# Note that for GRanges object, we need to specify seqnames, ranges and strand
# Others would be metacolumns


# You can glimpse of exported versions of this data with

export(mtar, "mtar01.bed")
cat(readLines('mtar01.bed',n=5), sep='\n')
#  chr12	8985196	  8985217	  .	0	-
#  chr7	  117095439	117095461	.	0	+
#  chr17	23750063	23750088	.	0	+
#  chr7	  27187934	27187957	.	0	-
#  chr17	43458622	43458643	.	0	-

mtar[1:5]
#       seqnames              ranges strand |         targets        mirname
#          <Rle>           <IRanges>  <Rle> |        <factor>       <factor>
#   [1]    chr12     8985197-8985217      - | ENST00000000412     hsa-miR-16
#   [2]     chr7 117095440-117095461      + | ENST00000003084 hsa-miR-509-3p
#   [3]    chr17   23750064-23750088      + | ENST00000003834    hsa-miR-612
#   [4]     chr7   27187935-27187957      - | ENST00000006015 hsa-miR-423-3p
#   [5]    chr17   43458623-43458643      - | ENST00000006101   hsa-miR-125b


# From the video


library(ERBS)
data(package='ERBS') # Useful command
# Data sets in package 'ERBS':
# GM12878                             
# HepG2  

f1 = dir(system.file("extdata",package="ERBS"), full=TRUE)
f1
# [1] "C:/Users/cheun/Documents/R/win-library/3.5/ERBS/extdata/ENCFF001VEH.narrowPeak"
# [2] "C:/Users/cheun/Documents/R/win-library/3.5/ERBS/extdata/ENCFF001VKD.narrowPeak"
# [3] "C:/Users/cheun/Documents/R/win-library/3.5/ERBS/extdata/sampleInfo.txt" 

readLines(f1[1], 4) # to read the first narrowPeak file
# [1] "chrX\t1509354\t1512462\t5\t0\t.\t157.92\t310\t32.000000\t1991"     "chrX\t26801421\t26802448\t6\t0\t.\t147.38\t310\t32.000000\t387"   
# [3] "chr19\t11694101\t11695359\t1\t0\t.\t99.71\t311.66\t32.000000\t861" "chr19\t4076892\t4079276\t4\t0\t.\t84.74\t310\t32.000000\t1508"

# The narrow peaks file shown is textual, with no headers


# Converting the narrowPeak to bedGraph format
# Recognizing the connection between the narrowPeak and bedGraph formats, we can import immediately to a GRanges.
# (this is for using others work in the past for annotations)

library(rtracklayer)
imp = import(f1[1], format="bedGraph")
imp
# GRanges object with 1873 ranges and 7 metadata columns:
#      seqnames            ranges strand |     score       NA.      NA.1      NA.2      NA.3      NA.4      NA.5
#         <Rle>         <IRanges>  <Rle> | <numeric> <integer> <logical> <numeric> <numeric> <numeric> <integer>
#  [1]     chrX   1509355-1512462      * |         5         0      <NA>    157.92       310        32      1991
# ...

export(imp, "demoex.bed")
cat(readLines('demoex.bed',n=5), sep='\n')
# chrX	1509354	  1512462	  .	5	.
# chrX	26801421	26802448	.	6	.
# chr19	11694101	11695359	.	1	.
# chr19	4076892	  4079276	  .	4	.
# chr3	53288567	53290767	.	9	.

# One could also convert it into many other formats

export(imp, 'demoex.gff3')
cat(readLines('demoex.gff3',n=5), sep='\n')
# ##gff-version 3
# ##source-version rtracklayer 1.40.3
# ##date 2018-09-05
# chrX	rtracklayer	sequence_feature	1509355	1512462	5	.	.	NA.=0;NA.2=157.92;NA.3=310;NA.4=32;NA.5=1991
# chrX	rtracklayer	sequence_feature	26801422	26802448	6	.	.	NA.=0;NA.2=147.38;NA.3=310;NA.4=32;NA.5=387


# So the basic message of this little section
# is that our data can become someone else's annotation.
# And someone else's experiments will become our annotation.


###


# AnnotationHub assessment
# and Others techniques in the intro video

# The AnnotationHub package can be used to obtain GRanges or 
# other suitably designed containers for institutionally curated annotation.


# Hub exploration via query()

library(BiocInstaller)

biocLite('AnnotationHub')
library(AnnotationHub)
ah = AnnotationHub()
mah = mcols(ah)
names(mah)
#  [1] "title"              "dataprovider"       "species"            "taxonomyid"         "genome"            
#  [6] "description"        "coordinate_1_based" "maintainer"         "rdatadateadded"     "preparerclass"     
# [11] "tags"               "rdataclass"         "rdatapath"          "sourceurl"          "sourcetype"  


# Looking at the dataprovider

sort(table(ah$dataprovider), decreasing = T)[1:10]
#                        BroadInstitute                               Ensembl 
#                                 18248                                 12003 
#                                  UCSC ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/ 
#                                 10337                                  1018 
#                              Haemcode                           Inparanoid8 
#                                   945                                   268 
#                                 Pazar                               Gencode 
#                                    91                                    28 
#                           UWashington                              Stanford 
#                                    25                                    24


# Query for data

query(ah,'glio')
# AnnotationHub with 145 records
# ...
# # retrieve records with, e.g., 'object[["AH22400"]]' 
# title                                                                   
# AH22400 | wgEncodeAwgDnaseDukeGlioblaUniPk.narrowPeak.gz                          
# AH23117 | wgEncodeAwgTfbsUtaGlioblaCtcfUniPk.narrowPeak.gz 
# ...


# checking on e.g. AH22400

mcols(query(ah,'glio'))[1,]
# Which shows the descriptions for the AH22400


# Retrieving record

AH22400 = query(ah,'glio')[['AH22400']]
AH22400
# GRanges object with 167960 ranges and 6 metadata columns:
#           seqnames                 ranges strand |        name     score
#              <Rle>              <IRanges>  <Rle> | <character> <numeric>
#       [1]     chr1       [564506, 564655]      * |        <NA>         0
#       [2]     chr1       [565506, 565655]      * |        <NA>         0
# ...

names(metadata(AH22400))
# [1] "AnnotationHubName" "File Name"         "Data Source"       "Provider"         
# [5] "Organism"          "Taxonomy ID"      


# getting the ID of the queried records

glio_ah = rownames(mcols(query(ah,'glio')))
glio_ah
#   [1] "AH22400" "AH23117" "AH23118" "AH24607" "AH24608" "AH24609" "AH24610" "AH25439"
#   [9] "AH25440" "AH25514" "AH25596" "AH29284" "AH29285" "AH29286" "AH29287" "AH29288"


# Other techniques in the assessments

sort(table(mah$species), decreasing=TRUE)[1:10]
#            Homo sapiens            Mus musculus Drosophila melanogaster              Bos taurus         Pan troglodytes 
#                   26032                    2416                     434                     328                     311 
#       Rattus norvegicus             Danio rerio           Gallus gallus   Monodelphis domestica             Felis catus 
#                     303                     297                     270                     258                     251

# We can see the number of entries devoted to H. sapiens and other species. 
# Queries about the content of the hub can be nested. Use

names(query(query(ah, "HepG2"), "CTCF"))
#  [1] "AH22249" "AH22531" "AH22693" "AH23134" "AH23186" "AH23305" "AH25111" "AH25112" "AH25461" "AH27539" "AH27540"
# [12] "AH27541" "AH27542"

# to find names of files that address binding of CTCF 
# (a transcription factor involved in chromatin structure regulation) to DNA from HepG2 cells.

# How many entries address CTCF binding in HepG2?:
length(names(query(query(ah, "HepG2"), "CTCF")))
# [1] 13


###


# OrgDb assessment
# Interactive tables for genomic annotation: assessment


# Counting genes in a cytoband

# Use org.Hs.eg.db to determine the number of different genes 
# in cytoband 17q21.1. 
# You would use the MAP field as a key, 
# so the code would involve select(org.Hs.eg.db, key="17q21.1", ... 
# and you must specify additional parameter settings to solve the problem.

# How many genes are present on 17q21.1?

library(BiocInstaller)
biocLite('org.Hs.eg.db')
library(org.Hs.eg.db)
# org - organism level; Hs - Homo Sapiens; eg - entre gene; SQLite based package

columns(org.Hs.eg.db)
#  [1] "ACCNUM"       "ALIAS"        "ENSEMBL"      "ENSEMBLPROT"  "ENSEMBLTRANS"
#  [6] "ENTREZID"     "ENZYME"       "EVIDENCE"     "EVIDENCEALL"  "GENENAME"    
# [11] "GO"           "GOALL"        "IPI"          "MAP"          "OMIM"        
# [16] "ONTOLOGY"     "ONTOLOGYALL"  "PATH"         "PFAM"         "PMID"        
# [21] "PROSITE"      "REFSEQ"       "SYMBOL"       "UCSCKG"       "UNIGENE"     
# [26] "UNIPROT"     

genes = select(org.Hs.eg.db, keys = '17q21.1', keytype = 'MAP', columns = 'GENENAME')
genes
#        MAP                                      GENENAME
# 1  17q21.1                   colony stimulating factor 3
# 2  17q21.1          proteasome 26S subunit, non-ATPase 3
# ...
dim(genes)
# [1] 22  2


# Gene Ontology categories in a cytoband

# Add the column type 'GO' to the select command 
# used to solve the previous problem.

# Select GO tags occurring among the five most common annotations 
# for genes on 17q21.1 (there may be more than one)

GO_gene = select(org.Hs.eg.db, keys = '17q21.1', keytype = 'MAP', columns = 'GO')
GO_gene
#         MAP         GO EVIDENCE ONTOLOGY
# 1   17q21.1 GO:0005125      IDA       MF
# 2   17q21.1 GO:0005125      NAS       MF
# 3   17q21.1 GO:0005130      IBA       MF
# ...
sort(table(GO_gene$GO),decreasing=T)[1:5]
# GO:0005829 GO:0000122 GO:0001222 GO:0005125 GO:0005576 
#          3          2          2          2          2 


# Traceable author statement annotation for a GO annotation

# How many of the GO annotations for ORMDL3 have TAS 
# (traceable author statement) as their evidence code?

ORMDL3 = select(org.Hs.eg.db, keys = 'ORMDL3', keytype = 'SYMBOL', columns = 'GO')
# Note that keytype = 'SYMBOL'

sort(table(ORMDL3$EVIDENCE),decreasing = T)
# TAS IDA IMP IBA IEA IPI 
#   4   2   2   1   1   1


# Enumerating genes where ER binds at TSS

# We'll start out by generating the addresses of all genes for Homo sapiens, 
# and then acquire the ER binding sites in liver cells.

library(Homo.sapiens)
g = genes(Homo.sapiens)
library(ERBS)
data(HepG2)

# Interpret the following computation:
kp = g[resize(g,1) %over% HepG2]

# kp is the set of genes whose transcription start site 
# lies in a HepG2 binding site.


# Creating and filtering an HTML5 report on gene annotation

# The DT package provides very nice interactive tabulation. 
# We'll generate a data.frame and publish it to the browser.

nn = names(kp)
m = select(Homo.sapiens, keys=nn, keytype="ENTREZID",
           columns=c("SYMBOL", "GENENAME", "TERM", "GO"))
# Note the keytype 'ENTREZID'

biocLite('DT')
library(DT)
datatable(m)

# datatable() will open a browser page with the table; use the search facility. 
# How many entries mention apoptosis?:

# Showing 1 to 10 of 31 entries (filtered from 1,084 total entries)


# Other techniques in the intro video


# Gene Ontology (GO) 

# is a widely used structured vocabulary that organizes terms relevant to the roles 
# of genes and gene products in biological processes, molecular functions, and
# cellular components.

library(BiocInstaller)
biocLite('GO.db')
library(GO.db)

# Looking at the keys in GO.db

keys(GO.db, keytype = 'TERM')[1:5]
# [1] "mitochondrion inheritance"                  
# [2] "mitochondrial genome maintenance"           
# [3] "reproduction"                               
# [4] "ribosome biogenesis"                        
# [5] "protein binding involved in protein folding"
GO_keys = keys(GO.db, keytype = 'TERM')

# Columns in GO.db
columns(GO.db)
# [1] "DEFINITION" "GOID"       "ONTOLOGY"   "TERM" 


# For example we would like to look at 'ribosome' 

ind = grep('ribosome', GO_keys)
GO_keys[ind]
#  [1] "ribosome biogenesis"                                                      
#  [2] "organellar ribosome"                                                      
#  [3] "structural constituent of ribosome" 
# ...

# And we are interested in 'ribosome biogenesis'

id = select(GO.db, keys = 'ribosome biogenesis', keytype = 'TERM', columns = 'GOID')
id
#                  TERM       GOID
# 1 ribosome biogenesis GO:0042254

# Knowing the GOID gives us a lot
# And we can search that in the org.Hs.eg.db

select(org.Hs.eg.db, keys = 'GO:0042254', keytype = 'GO', columns = 'SYMBOL') # 'GO' not 'GOID'
#            GO EVIDENCE ONTOLOGY  SYMBOL
# 1  GO:0042254      IBA       BP    GNL1
# 2  GO:0042254      IMP       BP     NVL
# 3  GO:0042254      IBA       BP   RPL34
# ...


###


# Assessment on KEGG


# Size of KEGG's human pathway catalog

# Use KEGGREST package (you need to have a live internet connection) 
# and obtain the vector returned by keggGet("hsa", "pathway")

# How many pathways are listed?

library(BiocInstaller)
biocLite('KEGG.db')
biocLite('KEGGREST')
# Note that the KEGGREST would be the updated one

library(KEGGREST)

# from forum:
# To obtain the number the KEGG's human pathway's, 
# run that command: keggList("pathway", "hsa"). 
# Count the number of pathway and remove 8. 
# (possible explanation: 8 new pathways have been added 
# since this assessment have been made)

hsa_pathway = keggList('pathway','hsa')
length(hsa_pathway)
# [1] 330


# Obtaining a pathway layout

# Continue with KEGGREST, and issue the commands

oo = keggGet("hsa00790", "image")

install.packages('png',dependencies = T)
library(png)
writePNG(oo, "im.png")

# Use a PNG viewer to examine the file im.png.
# What is the name of the pathway displayed?

# FOLATE BIOSYNTHESIS


# Other techniques from the intro video


# Finding BRCA2

library(org.Hs.eg.db)
columns(org.Hs.eg.db)

select(org.Hs.eg.db, keys = 'BRCA2', keytype = 'SYMBOL', columns = 'ENTREZID')
#   SYMBOL ENTREZID
# 1  BRCA2      675

library(KEGGREST)
brca2K = keggGet("hsa:675")
names(brca2K[[1]])
#  [1] "ENTRY"      "NAME"       "DEFINITION" "ORTHOLOGY"  "ORGANISM"   "PATHWAY"   
#  [7] "DISEASE"    "BRITE"      "POSITION"   "MOTIF"      "DBLINKS"    "STRUCTURE" 
# [13] "AASEQ"      "NTSEQ" 


# access to the columns

brca2K[[1]]$DISEASE
#                  H00019                  H00027                  H00031 
#     "Pancreatic cancer"        "Ovarian cancer"         "Breast cancer" 
#                  H00238                  H01554 
#        "Fanconi anemia" "Fallopian tube cancer" 


# Pathways in KEGG

hsa_pathway = keggList('pathway','hsa')

# For example I would like to search for fanconi anemia

grep('Fanconi',hsa_pathway)
# [1] 114
hsa_pathway[114]
#                                   path:hsa03460 
# "Fanconi anemia pathway - Homo sapiens (human)"

# then we can get the code of the pathway


# Getting the pathway from the KEGG

fanconi = keggGet('path:hsa03460')
names(fanconi[[1]])
#  [1] "ENTRY"       "NAME"        "DESCRIPTION" "CLASS"       "PATHWAY_MAP"
#  [6] "MODULE"      "DISEASE"     "DBLINKS"     "ORGANISM"    "GENE"       
# [11] "KO_PATHWAY"  "REFERENCE" 

fanconi_im = keggGet('path:hsa03460','image')
library(png)
writePNG(fanconi_im,'plot23.png')
# plot23


###


# Ontology lookup assessment


library(BiocInstaller)
biocLite('rols')
library(rols)
oo = Ontologies()

class(oo)
# [1] "Ontologies"
# attr(,"package")
# [1] "rols"

oo[[1]]
# Ontology: Agronomy Ontology (agro)  
# Ontology of agronomic practices, agronomic techniques, and agronomic
# variables used in agronomic experiments
# Loaded: 2018-05-15 Updated: 2018-10-02 Version: 2018-05-14 
# 1684 terms  713 properties  284 individuals


# Using rols to explore ontologies: a common term

# Load the rols library, and issue the command (with a live internet connection):
  
diab = OlsSearch("diabetes")

# What is the value of olsRows(allRows(diab))?

olsRows(allRows(diab))
# [1] 7478


# Continuing on, we can obtain a data.frame with all the lookup results.

fulld = olsSearch(allRows(diab))
adf = as(fulld, "data.frame")
names(adf)
#  [1] "id"                   "iri"                  "short_form"          
#  [4] "obo_id"               "label"                "description"         
#  [7] "ontology_name"        "ontology_prefix"      "type"                
# [10] "is_defining_ontology"

sort(table(adf$ontology_name), decreasing=TRUE)[1:10]
#   clo  ncit mondo   cco   efo  gexo  reto  rexo  ordo  doid 
#  5683   428   197   141   117   113   113   113   109    68 

# 'clo' stands for cell line ontology. 
# How many entries does it contribute on the term 'diabetes'?


# Searching an interactive table

# Continuing with the adf computed previously, use

library(DT)
datatable(adf)

# Search for 'oral glucose'. How many entries are found in the table?
# Showing 1 to 7 of 7 entries (filtered from 7,478 total entries)


###