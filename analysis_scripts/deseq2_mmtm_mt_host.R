library( "DESeq2" )
library(ggplot2)

#########################################################
# Read in data
#########################################################

#load the human transcriptomics results
mmtm_mt_host_comp <- read.csv("/Users/ichaudry/Documents/OneDrive - University of Maryland School of Medicine/research/ravel_lab/MMTM/analysis/02_mt_host_analysis/MMTM_host_featureCounts_0520204.csv")
rownames(mmtm_mt_host_comp) <- mmtm_mt_host_comp$Geneid
mmtm_mt_host_comp <- mmtm_mt_host_comp[,-1]

#filter samples with too few reads
min_reads <- 1000000
mmtm_mt_host_comp <- mmtm_mt_host_comp[,colSums(mmtm_mt_host_comp) >= min_reads]

#load in the cst assingments per sample
#sample_csts <- read.csv('/Users/ichaudry/Documents/OneDrive - University of Maryland School of Medicine/research/ravel_lab/MMTM/analysis/cst_analysis/mmtm_cst_assignments_taxa_corr.csv',sep=',',row.names = 1, header = TRUE)
#sample_csts <- read.csv('/Users/ichaudry/Documents/OneDrive - University of Maryland School of Medicine/research/ravel_lab/MMTM/analysis/01_mt_microbe_analysis/cst_assignments/mmtm_mt_cst_assignments_renamed.csv',sep=',',row.names = 1, header = TRUE)
sample_csts <- read.csv('/Users/ichaudry/Documents/OneDrive - University of Maryland School of Medicine/research/ravel_lab/MMTM/analysis/01_mt_microbe_analysis/cst_assignments/mmtm_mt_cst_assingments_06122024_renamed.csv',sep=',',row.names = 1, header = TRUE)

#load the metadata for each sample
mmtm_metadata <- read.csv("/Users/ichaudry/Documents/OneDrive - University of Maryland School of Medicine/research/ravel_lab/MMTM/analysis/mmtm_metadata.csv")
row.names(mmtm_metadata) <- mmtm_metadata$IGS.LABEL

#mmtm_metadata$IGS.LABEL <- sub("^", "MT_", mmtm_metadata$IGS.LABEL)
#mmtm_metadata$IGS.LABEL <- sub("$", "_ZYMO", mmtm_metadata$IGS.LABEL)

rownames(mmtm_metadata) <- mmtm_metadata$IGS.LABEL
mmtm_metadata <- mmtm_metadata[,-1]

#filter out the metadata dataframe to only include the samples that have gene data
samples = colnames(mmtm_mt_host_comp)
mmtm_metadata <- mmtm_metadata[row.names(mmtm_metadata) %in% samples, , drop=FALSE]

#filter out the control samples from the composition data by only keeping the samples that also have metadata associated with them
samples = rownames(mmtm_metadata)
mmtm_mt_host_comp <- mmtm_mt_host_comp[, colnames(mmtm_mt_host_comp) %in% samples, drop=FALSE]

#merge the CST information with the metadata
mmtm_metadata <- merge(sample_csts, mmtm_metadata, by='row.names' )
mmtm_metadata$CST <- as.factor(mmtm_metadata$CST)
mmtm_metadata$Timepoint <- as.factor(mmtm_metadata$Timepoint)




##############################################################
# Comparing expression within bodysite between CST-IV vs non-CST-IV
##############################################################

#pull out only cervical or vaginal samples
mmtm_metadata <- mmtm_metadata[mmtm_metadata$Cervix.Vagina == 'V',]
samples = mmtm_metadata$Row.names
mmtm_mt_host_comp <- mmtm_mt_host_comp[, colnames(mmtm_mt_host_comp) %in% samples, drop=FALSE]

#make another metadata variable telling if the sample is CSTIV or not
mmtm_metadata$isCSTIV <- (mmtm_metadata$CST == 'IV-A') | (mmtm_metadata$CST == 'IV-B') | (mmtm_metadata$CST == 'IV-C')
mmtm_metadata$isCSTIV <- as.factor(mmtm_metadata$isCSTIV)

#drop all genes that have 0 abundance accross all the samples
mmtm_mt_host_comp <- subset(mmtm_mt_host_comp, !rowSums(mmtm_mt_host_comp == 0) == ncol(mmtm_mt_host_comp))

mmtm_mt_host_comp <- as.matrix(mmtm_mt_host_comp)+1
#mmtm_metadata <- as.matrix(mmtm_metadata)


dds <- DESeqDataSetFromMatrix(countData=mmtm_mt_host_comp, 
                              colData=mmtm_metadata, 
                              design=~Timepoint + isCSTIV)

dds <- DESeq(dds)



res <- results(dds, contrast = c('isCSTIV', 'TRUE', 'FALSE'))
head(results(dds, tidy=TRUE)) #let's look at the results table


res <- res[order(res$padj),]
head(res)

write.csv(res, '/Users/ichaudry/Documents/OneDrive - University of Maryland School of Medicine/research/ravel_lab/MMTM/analysis/02_mt_host_analysis/deseq2/mmtm_mt_host_genes_deseq2_isCSTIV_vagina.csv')

#reset par
par(mfrow=c(1,1))
# Make a basic volcano plot
with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-6,6)))

# Add colored points: blue if padj<0.01, red if log2FC>1 and padj<0.05)
with(subset(res, padj<.01 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(res, padj<.01 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))

#########################################################


#########################################################
# Comparing expression between cervix and vagina within CST-IV and non CST-IV
#########################################################

#make another metadata variable telling if the sample is CSTIV or not
mmtm_metadata$isCSTIV <- (mmtm_metadata$CST == 'IV-A') | (mmtm_metadata$CST == 'IV-B')
mmtm_metadata$isCSTIV <- as.factor(mmtm_metadata$isCSTIV)


#pull out CSTIV or non-CSTIV samples
mmtm_metadata <- mmtm_metadata[mmtm_metadata$isCSTIV == FALSE,]
samples = mmtm_metadata$Row.names
mmtm_mt_host_comp <- mmtm_mt_host_comp[, colnames(mmtm_mt_host_comp) %in% samples, drop=FALSE]


#drop all genes that have 0 abundance accross all the samples
mmtm_mt_host_comp <- subset(mmtm_mt_host_comp, !rowSums(mmtm_mt_host_comp == 0) == ncol(mmtm_mt_host_comp))

mmtm_mt_host_comp <- as.matrix(mmtm_mt_host_comp)+1
#mmtm_metadata <- as.matrix(mmtm_metadata)


dds <- DESeqDataSetFromMatrix(countData=mmtm_mt_host_comp, 
                              colData=mmtm_metadata, 
                              design=~Timepoint + Cervix.Vagina)

dds <- DESeq(dds)



res <- results(dds, contrast = c("Cervix.Vagina", "C", "V"), tidy=TRUE)
head(results(dds, tidy=TRUE)) #let's look at the results table


res <- res[order(res$padj),]
head(res)

write.csv(res, '/Users/ichaudry/Documents/OneDrive - University of Maryland School of Medicine/research/ravel_lab/MMTM/analysis/02_mt_host_analysis/deseq2/mmtm_mt_host_genes_deseq2_CV_isNotCSTIV.csv')

#reset par
par(mfrow=c(1,1))
# Make a basic volcano plot
with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-3,3)))

# Add colored points: blue if padj<0.01, red if log2FC>1 and padj<0.05)
with(subset(res, padj<.05 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(res, padj<.05 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))


#########################################################

###############################################################
# Comparing expression within bodysite with a given CST overtime
###############################################################

#pull out only cervical or vaginal samples
mmtm_metadata <- mmtm_metadata[mmtm_metadata$Cervix.Vagina == 'V',]
samples = mmtm_metadata$Row.names
mmtm_mt_host_comp <- mmtm_mt_host_comp[, colnames(mmtm_mt_host_comp) %in% samples, drop=FALSE]

#pull out samples from a given CST
target_cst = c('IV-A', 'IV-B', 'IV-C')
#target_cst = c('I', 'II', 'III', 'V')
mmtm_metadata <- mmtm_metadata[mmtm_metadata$CST %in% target_cst,]
samples = mmtm_metadata$Row.names
mmtm_mt_host_comp <- mmtm_mt_host_comp[, colnames(mmtm_mt_host_comp) %in% samples, drop=FALSE]

#drop all genes that have 0 abundance accross all the samples
mmtm_mt_host_comp <- subset(mmtm_mt_host_comp, !rowSums(mmtm_mt_host_comp == 0) == ncol(mmtm_mt_host_comp))

mmtm_mt_host_comp <- as.matrix(mmtm_mt_host_comp)+1
#mmtm_metadata <- as.matrix(mmtm_metadata)


dds <- DESeqDataSetFromMatrix(countData=mmtm_mt_host_comp, 
                              colData=mmtm_metadata, 
                              design=~Timepoint)

dds <- DESeq(dds)



res <- results(dds, contrast = c('Timepoint', '1', '2'))
head(results(dds, tidy=TRUE)) #let's look at the results table


res <- res[order(res$padj),]
head(res)

write.csv(res, '/Users/ichaudry/Documents/OneDrive - University of Maryland School of Medicine/research/ravel_lab/MMTM/analysis/02_mt_host_analysis/deseq2/mmtm_mt_host_genes_deseq2_CSTIV_vagina_overtime.csv')

#reset par
par(mfrow=c(1,1))
# Make a basic volcano plot
with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-6,6)))

# Add colored points: blue if padj<0.01, red if log2FC>1 and padj<0.05)
with(subset(res, padj<.01 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(res, padj<.01 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))





