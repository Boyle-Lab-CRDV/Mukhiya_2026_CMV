library('edgeR')
library('BiocParallel')
library('tidyverse')
library('statmod')##issues loading
library(glmmSeq)
library(limma)
library(edgeR)
library(tidyverse)
library(kableExtra)
library(DESeq2)
library(AnnotationDbi)
library(org.Hs.eg.db)


meta_cmv= read.delim(file = "data/metadata_cmv_monos.txt", sep = "\t",
                 header=T, row.names=1)

meta_age= read.delim(file = "~/Documents/Mamba_Age_RNAseq/data/metadata_Age_monos.txt", sep = "\t",
                 header=T, row.names=1)

#original dim meta
dim(meta_age)
# 20 9
#Selecting only adults
meta_age1 = meta_age %>% filter(grepl(pattern = "Adult", x = age_bool))
dim(meta_age1)
#10 9

## joining datasets meta
colnames(meta_age1)
colnames(meta_cmv)

meta_age2 <- meta_age1[,c(1,5,6,7)]
meta_cmv1 <- meta_cmv[,c(1,4,8,9)]

meta_cmv_age <- rbind(meta_age2, meta_cmv1)
#change rownames
rownames(meta_cmv_age) <- meta_cmv_age$uid


##Counts datasets
counts_cmv <- read.csv("/Volumes/Macintosh HD/Users/zuleip/Documents/CMV_bulkRNAseq/CMV/data/Boyle_-_CMV_RNASeq_HTSeq.genes.RawCounts.csv", header=T, row.names = 1)

counts_age<- read.csv("~/Documents/Mamba_Age_RNAseq/data/Boyle_-_RNASeq_NextSeq550_HTSeq.genes.RawCounts.csv", header=T, row.names = 1)

dim(counts_age)
#57905    47

dim(counts_cmv)
#57905    52

##joining dataseta
counts_age_cmv <- cbind(counts_cmv, counts_age)
dim(counts_age_cmv)
#57905    99

##Removing children
counts_age_cmv1 <- counts_age_cmv[, -grep("Child", colnames(counts_age_cmv))]
dim(counts_age_cmv1)
#57905    79
##Removing Vd2
counts_age_cmv2 <- counts_age_cmv1[, -grep("Vd2", colnames(counts_age_cmv1))]
dim(counts_age_cmv2)
#57905    69
##Removing NK
counts_age_cmv3 <- counts_age_cmv2[, -grep("_NK_", colnames(counts_age_cmv2))]
dim(counts_age_cmv3)
# 57905    46
##Removing unpaired samples
counts_age_cmv4 <- counts_age_cmv3[, -grep("H22_M_NS", colnames(counts_age_cmv3))]
dim(counts_age_cmv4)
#57905    45
counts_age_cmv5 <- counts_age_cmv4[, -grep("H25_M_pRBC", colnames(counts_age_cmv4))]
dim(counts_age_cmv5)
#57905    44

##counts
counts_age_cmv5$ENSEMBL <- rownames(counts_age_cmv5)
counts_age_cmv6 <- counts_age_cmv5 %>% filter(grepl(pattern = "protein_coding", x = Gene.Biotype))
##Removing duplicates for now. 
counts_age_cmv6 <- counts_age_cmv6[!(duplicated(counts_age_cmv6$Associated.Gene.Name)), ]

##change rownames
rownames(counts_age_cmv6) <- counts_age_cmv6$Associated.Gene.Name
head(counts_age_cmv6)
counts_age_cmv7 <- counts_age_cmv6[,c(8:27, 35:44)]
##Checking quality
counts <- counts_age_cmv7
meta <- meta_cmv_age

##Create DGEList object
d0 <- DGEList(counts)

##Calculate normalization factors:
d0 <- calcNormFactors(d0)
d0

##Filtering low-expressed genes

dim(d0)
keep <- rowSums(d0[[1]]) >= 10
d <- d0[keep,]
dim(d)

##Multidimentional scaling (MDS) plot
#Visualizes the differences between the expression profiles 
#of different samples in two dimensions
##study of origin
group <- factor(meta$group) 
group

pdf(file="graphs/MD_monos_agecmv_pergroup.pdf")
points <- c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19) 
colors <- rep(c("blue","darkgreen","red", "purple", "orange","brown"),2)
plotMDS(cpm(d[[1]]), col=colors[group], pch=points[group], xlab = "leading FC dim1", ylab = "leading FC dim2") 
legend("topleft",legend=levels(group),pch=points,col=colors,ncol=2)
dev.off()



##Multidimentional scaling (MDS) plot
#Visualizes the differences between the expression profiles 
#of different samples in two dimensions
group <- factor(meta$state_bool) 
group


pdf(file="graphs/MD_monos_agecmv_perstate.pdf")
points <- c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19) 
colors <- rep(c("blue","darkgreen","red", "purple"),2)
plotMDS(cpm(d[[1]]), col=colors[group], pch=points[group], xlab = "leading FC dim1", ylab = "leading FC dim2") 
legend("topleft",legend=levels(group),pch=points,col=colors,ncol=2)
dev.off()



##Plot with samples ID
sam0 <- as.data.frame(colnames(counts))
sam0 %>% extract(col = `colnames(counts)`, into = "samples", regex = "(H\\d+)", remove = FALSE) -> sam0
sam0 %>% extract(col = `colnames(counts)`, into = "state", regex = "(...$)", remove = FALSE) -> sam0
sam0$state <- gsub("_NS", "unstim", sam0$state)
sam0$state <- gsub("RBC", "stim", sam0$state)
sam0$plotid = paste0(sam0$samples, "_", sam0$state)
head(sam0)

sam0 <- sam0[c(1:20),]

##Plot with samples ID
sam1 <- as.data.frame(colnames(counts))
sam1 %>% extract(col = `colnames(counts)`, into = "samples", regex = "(RDH\\d+)", remove = FALSE) -> sam1
sam1 %>% extract(col = `colnames(counts)`, into = "state", regex = "(Exviv|4hpRB\\w)", remove = FALSE) -> sam1
sam1$state <- gsub("Exviv", "unstim", sam1$state)
sam1$state <- gsub("4hpRBC", "stim", sam1$state)
sam1$plotid = paste0(sam1$samples, "_", sam1$state)
head(sam1)

sam1<- sam1[c(21:30),]

sam <- rbind(sam0, sam1)



plot_counts <- counts
colnames(plot_counts) <- sam$samples
sams <- factor(sam$samples)
samt <- factor(sam$state)

pdf(file="graphs/MDS_age_cmv_persample.pdf")
points <- c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19) 
colors <- rep(c("blue","darkgreen","red","purple","orange","brown","salmon","black","turquoise","darkgrey"),2)
plotMDS(cpm(d[[1]]), col=colors[samt], pch=points[sams], xlab = "leading FC dim1", ylab = "leading FC dim2") 
legend("topleft",legend=levels(sams),pch=points,col=colors,ncol=4)
dev.off()




