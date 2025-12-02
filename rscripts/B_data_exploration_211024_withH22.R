#Author: Z Pava
#.libPaths("/Users/zuleip/miniconda3/envs/rnaseq-age/lib/R/library")

library('edgeR')
library('BiocParallel')
library('tidyverse')
library('statmod')##issues loading
library(here)
############INPUT###########
##Analysis WITH H22

##Getting the HTSeq raw data counts

counts0 <- read.csv(here("input/Boyle_-_CMV_RNASeq_HTSeq.genes.RawCounts.csv"), header=T, row.names = 1)

samplesname <- as.character(colnames(counts0[, c(8:52)]))
meta1<- as.data.frame(counts0[1:45, c(8:length(counts0))], row.names = samplesname)
meta1$iud <- rownames(meta1)
meta <- as.data.frame(meta1[1:45, 46], row.names = samplesname)
colnames(meta)[1] = "uid"

#Extracting meta data from uid
meta
meta %>% extract(col = uid, into = "samples", regex = "(H\\d+)", remove = FALSE) -> meta
meta %>% extract(col = uid, into = "celltype", regex = "(_.+_)", remove = FALSE) -> meta
meta %>% extract(col = uid, into = "state", regex = "(...$)", remove = FALSE) -> meta

meta$celltype <- gsub("_M_", "monocytes", meta$celltype)
meta$celltype <- gsub("_NK_", "natkillers", meta$celltype)

meta$state <- gsub("_NS", "unstim", meta$state)
meta$state <- gsub("RBC", "stim", meta$state)

meta$cmv <- NA
meta$cmv[which(str_detect(meta$samples, "H10"))] <- "pos"
meta$cmv[which(str_detect(meta$samples, "H18"))] <- "pos"
meta$cmv[which(str_detect(meta$samples, "H25"))] <- "pos"
meta$cmv[which(str_detect(meta$samples, "H41"))] <- "pos"
meta$cmv[which(str_detect(meta$samples, "H48"))] <- "pos"
meta$cmv[which(str_detect(meta$samples, "H51"))] <- "pos"
meta$cmv[which(str_detect(meta$samples, "H22"))] <- "neg"
meta$cmv[which(str_detect(meta$samples, "H23"))] <- "neg"
meta$cmv[which(str_detect(meta$samples, "H40"))] <- "neg"
meta$cmv[which(str_detect(meta$samples, "H44"))] <- "neg"
meta$cmv[which(str_detect(meta$samples, "H49"))] <- "neg"
meta$cmv[which(str_detect(meta$samples, "H50"))] <- "neg"

meta$sex <- NA
meta$sex[which(str_detect(meta$samples, "H10"))] <- "M"
meta$sex[which(str_detect(meta$samples, "H18"))] <- "M"
meta$sex[which(str_detect(meta$samples, "H25"))] <- "F"
meta$sex[which(str_detect(meta$samples, "H41"))] <- "M"
meta$sex[which(str_detect(meta$samples, "H48"))] <- "F"
meta$sex[which(str_detect(meta$samples, "H51"))] <- "F"
meta$sex[which(str_detect(meta$samples, "H22"))] <- "F"
meta$sex[which(str_detect(meta$samples, "H23"))] <- "F"
meta$sex[which(str_detect(meta$samples, "H40"))] <- "F"
meta$sex[which(str_detect(meta$samples, "H44"))] <- "M"
meta$sex[which(str_detect(meta$samples, "H49"))] <- "M"
meta$sex[which(str_detect(meta$samples, "H50"))] <- "M"

meta$age_y <- NA
meta$age_y[which(str_detect(meta$samples, "H10"))] <- 35
meta$age_y[which(str_detect(meta$samples, "H18"))] <- 42
meta$age_y[which(str_detect(meta$samples, "H25"))] <- 45
meta$age_y[which(str_detect(meta$samples, "H41"))] <- 36
meta$age_y[which(str_detect(meta$samples, "H48"))] <- 26
meta$age_y[which(str_detect(meta$samples, "H51"))] <- 33
meta$age_y[which(str_detect(meta$samples, "H22"))] <- 29
meta$age_y[which(str_detect(meta$samples, "H23"))] <- 36
meta$age_y[which(str_detect(meta$samples, "H40"))] <- 66
meta$age_y[which(str_detect(meta$samples, "H44"))] <- 36
meta$age_y[which(str_detect(meta$samples, "H49"))] <- ""
meta$age_y[which(str_detect(meta$samples, "H50"))] <- ""

meta$group <- factor(paste0(meta$cmv, sep = ".", meta$state))
meta = meta %>% mutate(
  state_bool = case_when(state == "stim" ~ "1",
                         state == "unstim" ~ "0"),
  cmv_bool = case_when(cmv == "neg" ~ "0",
                       cmv == "pos" ~ "1"),
  cmv_state_bool = case_when(cmv == "pos" & state == "stim" ~ "1",
                             cmv == "neg" & state == "stim" ~ "3",
                             cmv == "pos" & state == "unstim" ~ "0",
                             cmv == "neg" & state == "unstim" ~ "2")
)

meta$state_bool <- factor(meta$state_bool, levels = c(0,1), labels=c("unstim","stim") )
meta$cmv_bool <- factor(meta$cmv_bool, levels = c(0,1), labels=c("neg","pos") )


#Subsetting per celltype.

meta_NK = meta %>% filter(grepl(pattern = "natkillers", x = celltype))  
meta_monos = meta %>% filter(grepl(pattern = "monocytes", x = celltype))
dim(meta_NK)
#23 11
dim(meta_monos)
#22 11
##counts
counts0$ENSEMBL <- rownames(counts0)
counts1 <- counts0 %>% filter(grepl(pattern = "protein_coding", x = Gene.Biotype))
##Removing duplicates for now. 
counts1 <- counts1[!(duplicated(counts1$Associated.Gene.Name)), ]

##change rownames
rownames(counts1) <- counts1$Associated.Gene.Name

##Subsetting per cell type and sorting. rownames(metadata) must be in the 
#same order as colnames(counts)

####REVIEW####
counts_NK <- counts1[, grep("NK", colnames(counts1))]
##Order meta rows to coincide the columns in the counts file.
meta_NK <- meta_NK[c("H10_NK_NS", "H10_NK_pRBC", "H18_NK_NS", "H18_NK_pRBC", "H22_NK_pRBC", "H23_NK_NS", "H23_NK_pRBC","H25_NK_NS",
                           "H25_NK_pRBC", "H40_NK_NS", "H40_NK_pRBC", "H41_NK_NS", "H41_NK_pRBC","H44_NK_NS", "H44_NK_pRBC",
                           "H48_NK_NS", "H48_NK_pRBC", "H49_NK_NS", "H49_NK_pRBC", "H50_NK_NS", "H50_NK_pRBC", "H51_NK_NS",
                           "H51_NK_pRBC"), ]

counts_monos <- counts1[, grep("M_", colnames(counts1))]
##Order meta rows to coincide the columns in the counts file.
meta_monos <- meta_monos[c("H10_M_NS", "H10_M_pRBC", "H18_M_NS", "H18_M_pRBC", "H22_M_NS", "H23_M_NS", "H23_M_pRBC",
                           "H25_M_pRBC", "H40_M_NS", "H40_M_pRBC", "H41_M_NS", "H41_M_pRBC","H44_M_NS", "H44_M_pRBC",
                           "H48_M_NS", "H48_M_pRBC", "H49_M_NS", "H49_M_pRBC", "H50_M_NS", "H50_M_pRBC", "H51_M_NS",
                           "H51_M_pRBC"), ]


###### NK data######

###Quality check
counts <- counts_NK
meta <- meta_NK

##Create DGEList object
d0 <- DGEList(counts)

##Calculate normalization factors:
#calcNormFactors doesn’t normalize the data, it just calculates 
#normalization factors for use downstream.
d0 <- calcNormFactors(d0)
d0

####Filtering low-expressed genes####
#Removing genes that are lowly expressed
#Here we perform a minimal pre-filtering to keep only rows that have 
#at least 10 reads total. Note that more strict filtering to increase 
#power is automatically applied via independent filtering on the mean of 
#normalized counts within the results function

dim(d0)
keep <- rowSums(d0[[1]]) >= 10
d <- d0[keep,]
dim(d)


#A mean-difference plot (MD-plot) is a plot of log-intensity ratios (differences) 
#versus log-intensity averages (means). Useful to detect outliers.

for (i in 1:20) {
  pdf(file= paste("output_withH22/graphs/MD_NK_cmv_",i,".pdf", sep = ""))
  plotMD(cpm(d[[1]], log=TRUE), column=i, xlab = "Average log-expression",
         ylab = "Expression log-ratio (this sample vs others)")
  abline(h=0,col="red",lty=2,lwd=2)
  dev.off()
}

group <- factor(meta$group)
group

##Multidimentional scaling (MDS) plot
#Visualizes the differences between the expression profiles 
#of different samples in two dimensions

pdf(file="output_withH22/graphs/MDS_NK_cmv_pergroup.pdf")
points <- c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19)
colors <- rep(c("blue","darkgreen","red", "purple"),2)
plotMDS(d[[1]], col=colors[group], pch=points[group], xlab = "leading FC dim1", ylab = "leading FC dim2")
legend("topleft",legend=levels(group),pch=points,col=colors,ncol=2)
dev.off()

##Plot with samples ID
sam <- as.data.frame(colnames(counts))
sam %>% extract(col = `colnames(counts)`, into = "samples", regex = "(H\\d+)", remove = FALSE) -> sam
sam %>% extract(col = `colnames(counts)`, into = "state", regex = "(...$)", remove = FALSE) -> sam
sam$state <- gsub("_NS", "unstim", sam$state)
sam$state <- gsub("RBC", "stim", sam$state)
sam$plotid = paste0(sam$samples, "_", sam$state)
head(sam)

plot_counts <- counts
colnames(plot_counts) <- sam$plotid
sams <- factor(sam$plotid)
# 
pdf(file=here("output_withH22/graphs/MDS_NK_cmv_persample.pdf"))
points <- c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19)
colors <- rep(c("blue","darkgreen","red","purple","orange","brown","salmon","black","turquoise","darkgrey"),2)
plotMDS(d[[1]], col=colors[sams], pch=points[sams], xlab = "leading FC dim1", ylab = "leading FC dim2")
legend("topleft",legend=levels(sams),pch=points,col=colors,ncol=4)
dev.off()

##Saving counts and metadata
d[[1]] -> counts_NK
meta -> meta_NK

write.table(counts_NK, file = here("output_withH22/data/counts_cmv_NK_withH22_oct24.txt"), sep = "\t",
            row.names = TRUE)
write.table(meta_NK, file = here("output_withH22/data/metadata_cmv_NK_withH22_oct24.txt"), sep = "\t",
            row.names = TRUE)


# ######Monocyte data######
# 
###Quality check
counts <- counts_monos
meta <- meta_monos

##Create DGEList object
d0 <- DGEList(counts)

##Calculate normalization factors:
#calcNormFactors doesn’t normalize the data, it just calculates 
#normalization factors for use downstream.
d0 <- calcNormFactors(d0)
d0

####Filtering low-expressed genes####
#Removing genes that are lowly expressed
#Here we perform a minimal pre-filtering to keep only rows that have 
#at least 10 reads total. Note that more strict filtering to increase 
#power is automatically applied via independent filtering on the mean of 
#normalized counts within the results function

dim(d0)
# 20268    22
keep <- rowSums(d0[[1]]) >= 10
d <- d0[keep,]
dim(d)
# 15561    22

#A mean-difference plot (MD-plot) is a plot of log-intensity ratios (differences) 
#versus log-intensity averages (means). Useful to detect outliers.

for (i in 1:22) {
  pdf(file= paste("output_withH22/graphs/MD_monos_cmv_",i,".pdf", sep = ""))
  plotMD(cpm(d[[1]], log=TRUE), column=i, xlab = "Average log-expression",
         ylab = "Expression log-ratio (this sample vs others)")
  abline(h=0,col="red",lty=2,lwd=2)
  dev.off()
}

group <- factor(meta$group)
group

##Multidimentional scaling (MDS) plot
#Visualizes the differences between the expression profiles 
#of different samples in two dimensions

pdf(file="output_withH22/graphs/MDS_monos_cmv_pergroup.pdf")
points <- c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19)
colors <- rep(c("blue","darkgreen","red", "purple"),2)
plotMDS(d[[1]], col=colors[group], pch=points[group], xlab = "leading FC dim1", ylab = "leading FC dim2")
legend("topleft",legend=levels(group),pch=points,col=colors,ncol=2)
dev.off()

##Plot with samples ID
sam <- as.data.frame(colnames(counts))
sam %>% extract(col = `colnames(counts)`, into = "samples", regex = "(H\\d+)", remove = FALSE) -> sam
sam %>% extract(col = `colnames(counts)`, into = "state", regex = "(...$)", remove = FALSE) -> sam
sam$state <- gsub("_NS", "unstim", sam$state)
sam$state <- gsub("RBC", "stim", sam$state)
sam$plotid = paste0(sam$samples, "_", sam$state)
head(sam)

plot_counts <- counts
colnames(plot_counts) <- sam$plotid
sams <- factor(sam$plotid)
# 
pdf(file=here("output_withH22/graphs/MDS_monos_cmv_persamplewleg.pdf"))
points <- c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19)
colors <- rep(c("blue","darkgreen","red","purple","orange","brown","salmon","black","turquoise","darkgrey"),2)
plotMDS(d[[1]], col=colors[sams], pch=points[sams], xlab = "leading FC dim1", ylab = "leading FC dim2")
legend("topleft",legend=levels(sams),pch=points,col=colors,ncol=4)
dev.off()

##Saving counts and metadata
d[[1]] -> counts_monos
meta -> meta_mono

write.table(counts_monos, file = here("output_withH22/data/counts_cmv_monos_withH22_oct24.txt"), sep = "\t",
            row.names = TRUE)
write.table(meta_mono, file = here("output_withH22/data/metadata_cmv_monos_withH22_oct24.txt"), sep = "\t",
            row.names = TRUE)
##done 240725
writeLines(capture.output(sessionInfo()), "session_info_240725.txt")
