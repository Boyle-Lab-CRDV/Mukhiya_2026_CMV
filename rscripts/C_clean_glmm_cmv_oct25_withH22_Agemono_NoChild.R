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

#### Age Monocytes analysis ####

#### INPUTS ####

## The monocytes experiment for the cmv project failed. 
## In this analysis, we are aiming to investigate the effect of cmv
## in the age monocyte dataset. Previously we had demonstrated that
## age has an effect modifier in the response to malaria.
## The current model will take that into account


#Reading clean expected counts data see script "data_exploration.R" to
#reproduce these files.
counts= read.delim(file = "/working_groups/boylelab/shared/ZulyPava/cmv_project_R4.2/output_withH22/data/counts_Age_monos_cmv.txt", sep = "\t",
                   header=T, row.names=1)
meta= read.delim(file = "/working_groups/boylelab/shared/ZulyPava/cmv_project_R4.2/output_withH22/data/metadata_Age_monos_cmv.txt", sep = "\t",
                 header=T, row.names=1)

head(counts)
head(meta)
# 
dim(counts)
#15274    20

dim(meta)
#20 11
table(colnames(counts) %in% rownames(meta))

colnames(meta)[10]<-"cmv_bool"
colnames(meta)[11]<-"cmv"
## Data formatting
##removing children from the dataset
meta <- meta%>%
  filter(age=="Adult")
counts <- counts[, rownames(meta)]

###MODEL FITTING#####
#Using negative binomial models requires gene dispersion estimates to be made.
#Calculating dispersion with DESeq.
meta$state_bool <- factor(meta$state_bool, levels = c("unstim","stim"), labels=c("unstim","stim") )
colnames(meta)[10]="cmv_bool"
colnames(meta)[11]="cmv"

# meta = meta%>% mutate(cmv_bool=case_when( cmv=="1"~"pos",
#                                           cmv=="0"~"neg",
#                                           TRUE~NA_character_)) 
#meta$age_bool <- factor(meta$age_bool, levels = c("Adult", "Child"), labels=c("Adult", "Child") )
meta$cmv_bool <- factor(meta$cmv_bool, levels = c("Negative","Positive"), labels=c("Negative","Positive") )
meta$cmv <- factor(meta$cmv, levels = c("0","1"), labels=c("0","1") )
##matching order of colnames and rownames
counts <- counts[,rownames(meta)]

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

pdf(file = here("output_withH22/graphs/monocytes_ageNochild_cmv_Dispersal estimates_DESeq2.pdf"))
plotDispEsts(dds, main="monocytes Dispersal estimates_DESeq2")
dev.off()

sum(res$padj < 0.05, na.rm=TRUE)

pdf(file = here("output_withH22/graphs/monocytes_ageNochild_cmv_MA Plot_DESeq2.pdf"))
plotMA(res, ylim=c(-2,2), main="monocytes_age_cmv_MA Plot_DESeq2")
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
glmm_results3 <- glmmSeq(~state_bool*cmv_bool + (1 | samples),
                        countdata = counts,
                        metadata = meta,
                        dispersion = dispersions,
                        progress=TRUE,
                        cores = 32)
#Errors in 4 gene(s): NAALAD2, ZNF831, CDHR3, GYLTL1B
glmm_results3 <- glmmQvals(glmm_results3, cutoff = 0.05, verbose = TRUE) 
saveRDS(glmm_results3, "/working_groups/boylelab/shared/ZulyPava/cmv_project_R4.2/output_withH22/data/glmm_results3.rds")

# state_bool
# ----------
#   Not Significant     Significant 
# 7413            7851 
# 
# cmv_bool
# --------
#   Not Significant     Significant 
# 13487            1777 
# 
# state_bool:cmv_bool
# -------------------
#   Not Significant     Significant 
# 13380            1884 

#### Adjust to the actual data from here!!! ####
#### RESULTS ####

##Getting data for plots
stats = data.frame(glmm_results3@stats)
predict = data.frame(glmm_results3@predict)
stats$gene.id <- rownames(stats)
predict$gene.id <- rownames(predict)

##qval similar to p.adjust
plotdata0 <- left_join(predict, stats, by = "gene.id")
plotdata0$ENSEMBL <- plotdata0$gene.id

##Changing the ensembl ids for gene names
keys(org.Hs.eg.db, keytype="ENSEMBL")[1:10]

anno <- AnnotationDbi::select(org.Hs.eg.db,keys=plotdata0$ENSEMBL,
                              columns=c("ENSEMBL","SYMBOL","GENENAME"),
                              keytype="ENSEMBL")
# Have a look at the annotation
head(anno)
dim(anno)
#15358     3

dup_ids <- anno$ENSEMBL[duplicated(anno$ENSEMBL)]

filter(anno, ENSEMBL %in% dup_ids) %>% 
  arrange(ENSEMBL) %>% head

anno <- AnnotationDbi::select(org.Hs.eg.db,keys=plotdata0$ENSEMBL,
                              columns=c("ENSEMBL","SYMBOL","GENENAME"),
                              keytype="ENSEMBL") %>% filter(!duplicated(ENSEMBL))

dim(anno)
# 15264     3
##Calculating FC and subsetting for plotting
plotdata <- left_join(plotdata0, anno,by="ENSEMBL")
colnames(plotdata)
##Fixed the direction of the fold change 16-July-2025
# [1] "y_unstim_Negative"                      "y_stim_Negative"                        "y_unstim_Positive"                     
# [4] "y_stim_Positive"  

plotdata1 = plotdata %>% mutate(
  #State_cmvnFC=(y_stim_Negative_Adult-y_unstim_Negative_Adult)
  State_cmvnFC_Adults = log2(plotdata[, 2]+1) - log2(plotdata[, 1]+1),
  #State_cmvpFC=(y_stim_Positive_Adult-y_unstim_Positive_Adult)
  State_cmvpFC_Adults = log2(plotdata[, 4]+1) - log2(plotdata[, 3]+1),
  #CMV_unstimFC=(y_unstim_Positive_Adult-y_unstim_Negative_Adult)
  CMV_unstimFC_Adults = log2(plotdata[, 3]+1) - log2(plotdata[, 1]+1),
  #CMV_stimFC=(y_stim_Positive_Adult-y_stim_Negative_Adult)
  CMV_stimFC_Adults = log2(plotdata[, 4]+1) - log2(plotdata[, 2]+1)) 

colnames(plotdata1)

plotdata2 <- plotdata1[,c(38:44,32:37)]
##Coding for variable SigG
## cmv = Sig for cmv but NS for state or state:cmv
## state = Sig for state but NS for cmv or state:cmv
## state:cmv = Sig for state:cmv with interaction
## state&cmv = Sig for state and cmv but without interaction.
## ns = NS for cmv, state or state:cmv

alpha <- 0.05

plotdata2 <- plotdata2 %>%
  mutate(
    # helper flags
    s  = qvals.state_bool               <= alpha,   # state main
    c  = qvals.cmv_bool                 <= alpha,   # cmv main
    sc = qvals.state_bool.cmv_bool      <= alpha,   # state:cmv interaction

    SigG = case_when(
      # --- interactions take precedence ---
      sc                             ~ "state&cmv",               # state:cmv interaction sig (regardless of mains)
      # --- main-effect combinations (no interactions) ---
      s & c                          ~ "state:cmv",               # two mains sig
      s                              ~ "state_only",              # single main sig
      c                              ~ "cmv_only",

      # --- none significant ---
      !(s | c | sc)         ~ "ns",
      
      # --- safety net for odd NA combos ---
      TRUE                           ~ "unclassified"
    )
  )

table(plotdata2$SigG, useNA = "ifany")

# cmv_only         ns state_only  state:cmv  state&cmv 
# 162       6746       5896        576       1884 
##Check SigG coding is correct
tapply(plotdata2$qvals.state_bool.cmv_bool, plotdata2$SigG, summary)
tapply(plotdata2$qvals.state_bool, plotdata2$SigG, summary)
tapply(plotdata2$qvals.cmv_bool, plotdata2$SigG, summary)


##Getting significant gene lists
plotdata3 <- subset(plotdata2, SigG != "ns")

##Saving results
write.table(plotdata1, file = "output_withH22/data/glmmseq_cmv_AgeMonocytes_NoChild_plot1_oct25.txt", sep = "\t",
            row.names = TRUE)
write.table(plotdata2, file = "output_withH22/data/glmmseq_cmv_AgeMonocytes_NoChild_plot2_oct25.txt", sep = "\t",
            row.names = TRUE)
write.table(plotdata3, file = "output_withH22/data/glmmseq_cmv_AgeMonocytes_NoChild_plot3_oct25.txt", sep = "\t",
            row.names = TRUE)

write.csv(plotdata1, file = "output_withH22/data/glmmseq_cmv_AgeMonocytes_NoChild_plot1_oct25.csv")
write.csv(plotdata2, file = "output_withH22/data/glmmseq_cmv_AgeMonocytes_NoChild_plot2_oct25.csv")
write.csv(plotdata3, file = "output_withH22/data/glmmseq_cmv_AgeMonocytes_NoChild_plot3_oct25.csv")

##done 240725
writeLines(capture.output(sessionInfo()), "session_info_age_mon_101125.txt")


