#' ---
#' title: "DE analysis using glmm of NK and monocytes from cmv+/- people stimulated with P. falciparum"
#' Date: 20 Dec 2021
#' output: html_document
#' --- 
#'
#' Notes:glmmSeq allows model to include random effects, like patient ID.
#' 
#' 


library(glmmSeq)
library(limma)
library(edgeR)
library(tidyverse)
library(kableExtra)
library(DESeq2)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(here)

#### NK analysis ####

#### INPUTS ####

#Reading clean expected counts data see script "data_exploration.R" to
#reproduce these files.

counts = read.delim(file = here("output_withH22/data/counts_cmv_NK_withH22_oct24.txt"), sep = "\t",
                        header=T, row.names=1)
head(counts)
meta = read.delim(file = here("output_withH22/data/metadata_cmv_NK_withH22_oct24.txt"), sep = "\t",
                      header=T, row.names=1)
head(meta)
# 
dim(counts)
#15206    22
# 15263    23 #with H22
dim(meta)
#22 11
# 23 11 with H22


###MODEL FITTING#####
#Using negative binomial models requires gene dispersion estimates to be made.
#Calculating dispersion with DESeq.
meta$state_bool <- factor(meta$state_bool, levels = c("unstim","stim"), labels=c("unstim","stim") )
meta$cmv_bool <- factor(meta$cmv_bool, levels = c("neg","pos"), labels=c("neg","pos") )


dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = meta,
                              design = ~state_bool*cmv_bool)

dds <- DESeq(dds)

dispersions <- setNames(dispersions(dds), rownames(counts))

#DESeq2 dispersion estimates are inversely related to the mean and directly 
#related to variance. Based on this relationship, the dispersion is higher 
#for small mean counts and lower for large mean counts.
res <- results(dds)
res

res <- res[order(res$padj),]
res

pdf(file = here("output_withH22/graphs/NK Dispersal estimates_DESeq2.pdf"))
plotDispEsts(dds, main="NK Dispersal estimates_DESeq2")
dev.off()

sum(res$padj < 0.05, na.rm=TRUE)

pdf(file = here("output_withH22/graphs/NK MA Plot_DESeq2.pdf"))
plotMA(res, ylim=c(-2,2), main="NK MA Plot_DESeq2")
dev.off()

#An MA-plot is a scatter plot of log2 fold changes (on the y-axis) versus the mean of normalized counts (on the x-axis)
## The MA-plot shows the log2 fold changes attributable to a given variable over the mean of normalized counts 
## for all the samples in the DESeqDataSet. Points will be colored red if the adjusted p value is less than 0.1. 
## Points which fall out of the window are plotted as open triangles pointing either up or down.


##Estimating size factors 
sizeFactors <- estimateSizeFactorsForMatrix(counts)

###MODEL###
set.seed(12345)
##Running the model state accounting for age
glmm_results <- glmmSeq(~state_bool*cmv_bool + (1 | samples),
                        countdata = counts,
                        metadata = meta,
                        dispersion = dispersions,
                        progress=TRUE,
                        cores = 32)
#Errors in 1 gene(s): RNF17
glmm_results <- glmmQvals(glmm_results, cutoff = 0.05, verbose = TRUE) 



##With H22

# state_bool
# ----------
#   Not Significant     Significant 
# 5695            9567 
# 
# cmv_bool
# --------
#   Not Significant     Significant 
# 14998             264 
# 
# state_bool:cmv_bool
# -------------------
#   Not Significant     Significant 
# 15151             111 

#### RESULTS ####
##Getting data for plots
stats = data.frame(glmm_results@stats)
predict = data.frame(glmm_results@predict)
stats$gene.id <- rownames(stats)
predict$gene.id <- rownames(predict)

##qval similar to p.adjust

plotdata0 <- left_join(predict, stats, by = "gene.id")
plotdata0$SYMBOL <- plotdata0$gene.id

##Changing the ensembl ids for gene names
keys(org.Hs.eg.db, keytype="SYMBOL")[1:10]

anno <- AnnotationDbi::select(org.Hs.eg.db,keys=plotdata0$SYMBOL,
                              columns=c("ENSEMBL","SYMBOL","GENENAME"),
                              keytype="SYMBOL")
# Have a look at the annotation
head(anno)
dim(anno)
#17213     3

dup_ids <- anno$ENSEMBL[duplicated(anno$ENSEMBL)]

filter(anno, ENSEMBL %in% dup_ids) %>% 
  arrange(ENSEMBL) %>% head

anno <- AnnotationDbi::select(org.Hs.eg.db,keys=plotdata0$SYMBOL,
                               columns=c("ENSEMBL","SYMBOL","GENENAME"),
                               keytype="SYMBOL") %>% filter(!duplicated(ENSEMBL))

dim(anno)
#15734      3
##Calculating FC and subsetting for plotting
plotdata <- left_join(plotdata0, anno,by="SYMBOL")
head(plotdata)
##Fixed the direction of the fold change 16-July-2025
# [1] "y_unstim_neg"                    [2]  "y_stim_neg"                       
# [3] "y_unstim_pos"                    [4]  "y_stim_pos"

plotdata1 = plotdata %>% mutate(
  #State_cmvnFC=(y_stim_neg-y_unstim_neg)
  State_cmvnFC = log2(plotdata[, 2]+1) - log2(plotdata[, 1]+1),
  #State_cmvpFC=(y_stim_pos-y_unstim_pos)
  State_cmvpFC = log2(plotdata[, 4]+1) - log2(plotdata[, 3]+1),
  #CMV_unstimFC=(y_unstim_pos-y_unstim_neg)
  CMV_unstimFC = log2(plotdata[, 3]+1) - log2(plotdata[, 1]+1),
  #CMV_stimFC=(y_stim_pos-y_stim_neg)
  CMV_stimFC = log2(plotdata[, 4]+1) - log2(plotdata[, 2]+1)) 

colnames(plotdata1)

plotdata2 <- plotdata1[,c(38:44,35:37)]
##Coding for variable SigG
## cmv = Sig for cmv but NS for state or state:cmv
## state = Sig for state but NS for cmv or state:cmv
## state:cmv = Sig for state:cmv with interaction
## state&cmv = Sig for state and cmv but without interaction.
## ns = NS for cmv, state or state:cmv
plotdata2$SigG[plotdata2$qvals.state_bool>= 0.05 & plotdata2$qvals.state_bool.cmv_bool>=0.05 & plotdata2$qvals.cmv_bool>=0.05 ] <- "ns"
plotdata2$SigG[plotdata2$qvals.state_bool>= 0.05 & plotdata2$qvals.state_bool.cmv_bool>=0.05 & plotdata2$qvals.cmv_bool<=0.05 ] <- "cmv"
plotdata2$SigG[plotdata2$qvals.state_bool<= 0.05 & plotdata2$qvals.state_bool.cmv_bool>=0.05 & plotdata2$qvals.cmv_bool>=0.05 ] <- "state"
plotdata2$SigG[plotdata2$qvals.state_bool<= 0.05 & plotdata2$qvals.state_bool.cmv_bool>=0.05 & plotdata2$qvals.cmv_bool<=0.05 ] <- "state&cmv"
plotdata2$SigG[plotdata2$qvals.state_bool.cmv_bool<0.05] <- "state:cmv"
table(plotdata2$SigG, useNA = "ifany")

##Without H22
#cmv        ns     state state:cmv state&cmv 
#37      5469      9440       138       120

##With H22
# cmv        ns     state state:cmv state&cmv 
# 86      6412     10386       122       175 

tapply(plotdata2$qvals.state_bool.cmv_bool, plotdata2$SigG, summary)

# $`state:cmv`
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.000000 0.000000 0.001618 0.010958 0.023802 0.046155 
# 
# $`state&cmv`
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.09638 0.48067 0.69468 0.64024 0.82957 0.87883

tapply(plotdata2$qvals.state_bool, plotdata2$SigG, summary)

# $state
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 0.000e+00 0.000e+00 4.006e-05 5.559e-03 4.390e-03 4.998e-02 
 
tapply(plotdata2$qvals.cmv_bool, plotdata2$SigG, summary)
# $cmv
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.000000 0.005468 0.021778 0.020584 0.029948 0.048089  

##Coding for variable SigG
plotdata2$SigGbyCMV[plotdata2$SigG == "state:cmv" & plotdata2$State_cmvnFC > plotdata2$State_cmvpFC] <- "blue_cmvneg_SG_state:cmv"
plotdata2$SigGbyCMV[plotdata2$SigG == "state:cmv" & plotdata2$State_cmvnFC < plotdata2$State_cmvpFC] <- "yellow_cmvpos_SG_state:cmv"
table(plotdata2$SigGbyCMV, useNA ="ifany")


# blue_cmvneg_SG_state:cmv yellow_cmvpos_SG_state:cmv                       <NA> 
#   111                         11                      17059

##Getting significant gene lists
plotdata3 <- subset(plotdata2, SigG != "ns")

##Saving results
write.table(plotdata1, file = "output_withH22/data/glmmseq_cmv_NK_output_withH22_oct24.txt", sep = "\t",
            row.names = TRUE)
write.table(plotdata2, file = "output_withH22/data/glmmseq_cmv_NK_output_summary_withH22_oct24.txt", sep = "\t",
            row.names = TRUE)
write.table(plotdata3, file = "output_withH22/data/glmmseq_cmv_NK_output_SiGeneList_withH22_oct24.txt", sep = "\t",
            row.names = TRUE)

write.csv(plotdata1, file = "output_withH22/data/glmmseq_cmv_NK_output_withH22_oct24.csv")
write.csv(plotdata2, file = "output_withH22/data/glmmseq_cmv_NK_output_summary_withH22_oct24.csv")
write.csv(plotdata3, file = "output_withH22/data/glmmseq_cmv_NK_output_SiGeneList_withH22_oct24.csv")


