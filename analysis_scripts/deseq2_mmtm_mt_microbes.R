
#############################################################
#imports
#############################################################

library( "DESeq2" )
library(ggplot2)
library('statmod')
library('dplyr')

#########################################################
#  Compare the relative kegg expression
#  between the cervical and vaginal samples using Deseq2
# (looking across all CSTs) but a specific timepoint
#########################################################

#read in kegg counts per sample
#mmtm_mt_comp <- read.csv("/Users/ichaudry/Library/CloudStorage/OneDrive-UniversityofMarylandSchoolofMedicine/research/ravel_lab/MMTM/analysis/01_mt_microbe_analysis/virgo2/taxa_filt_mmtm_mt_virgo2_kegg_counts_05132024.csv", header=TRUE, row.names=1)
mmtm_mt_comp <- read.csv("/Users/ichaudry/Library/CloudStorage/OneDrive-UniversityofMarylandSchoolofMedicine/research/ravel_lab/MMTM/analysis/01_mt_microbe_analysis/virgo2/taxa_filt_mmtm_mt_virgo2_kegg_counts_06122024_junk-genes-removed.csv", header=TRUE, row.names=1)

#filter samples with too few
min_reads <- 1000000
mmtm_mt_comp <- mmtm_mt_comp[,colSums(mmtm_mt_comp) >= min_reads]

#drop control columns
mmtm_mt_comp <- select(mmtm_mt_comp, starts_with("MT"))

#remove the MT_ and _ZYMO suffix from the column names
col_names <- colnames(mmtm_mt_comp)
col_names <- gsub("^MT_([A-Za-z0-9]+)_ZYMO$", "\\1", col_names)
colnames(mmtm_mt_comp) <- col_names

#load in the cst assingments per sample
#sample_csts <- read.csv('/Users/ichaudry/Documents/OneDrive - University of Maryland School of Medicine/research/ravel_lab/MMTM/analysis/cst_analysis/mmtm_cst_assignments_taxa_corr.csv',sep=',',row.names = 1, header = TRUE)
sample_csts <- read.csv('/Users/ichaudry/Documents/OneDrive - University of Maryland School of Medicine/research/ravel_lab/MMTM/analysis/01_mt_microbe_analysis/cst_assignments/mmtm_mt_cst_assingments_06122024_renamed.csv',sep=',',row.names = 1, header = TRUE)


#load the metadata for each sample
mmtm_metadata <- read.csv("/Users/ichaudry/Documents/OneDrive - University of Maryland School of Medicine/research/ravel_lab/MMTM/analysis/mmtm_metadata.csv", sep=',', row.names = 1)

#filter out the metadata dataframe to only include the samples that have composition data
samples = colnames(mmtm_mt_comp)
mmtm_metadata <- mmtm_metadata[row.names(mmtm_metadata) %in% samples, , drop=FALSE]
mmtm_metadata$Cervix.Vagina <- as.factor(mmtm_metadata$Cervix.Vagina)
mmtm_metadata$Timepoint <- as.factor(mmtm_metadata$Timepoint)
mmtm_metadata$Ext.Participant.ID <- as.factor(mmtm_metadata$Ext.Participant.ID)

#merge the CST information with the metadata
mmtm_metadata <- merge(sample_csts, mmtm_metadata, by='row.names' )
rownames(mmtm_metadata) <- mmtm_metadata$Row.names
mmtm_metadata$CST <- as.factor(mmtm_metadata$CST)

#filter out the control samples from the composition data by only keeping the samples that also have metadata associated with them
samples = rownames(mmtm_metadata)
mmtm_mt_comp <- mmtm_mt_comp[, colnames(mmtm_mt_comp) %in% samples, drop=FALSE]

#drop keggs with low abundance
mmtm_mt_comp <- mmtm_mt_comp[rowSums(mmtm_mt_comp[]) > 50,]

target_time = c('2')
mmtm_metadata <- mmtm_metadata[mmtm_metadata$Timepoint %in% target_time,]
samples = rownames(mmtm_metadata)
mmtm_mt_comp <- mmtm_mt_comp[, colnames(mmtm_mt_comp) %in% samples, drop=FALSE]

mmtm_mt_comp <- as.matrix(mmtm_mt_comp)+1
#mmtm_metadata <- as.matrix(mmtm_metadata)

dds <- DESeqDataSetFromMatrix(countData=mmtm_mt_comp, 
                              colData=mmtm_metadata, 
                              design=~CST + Cervix.Vagina)

dds <- DESeq(dds)



res <- results(dds, contrast = c("Cervix.Vagina", "C", "V"), tidy = TRUE)
head(res) #let's look at the results table


res <- res[order(res$padj),]
head(res)

write.csv(res, '/Users/ichaudry/Library/CloudStorage/OneDrive-UniversityofMarylandSchoolofMedicine/research/ravel_lab/MMTM/analysis/01_mt_microbe_analysis/deseq2/mmtm_mt_diff_micro_kegg_cervix-vagina_timepoint2_deseq2.csv', row.names = FALSE)

#reset par
par(mfrow=c(1,1))
# Make a basic volcano plot
with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-5,5)))

# Add colored points: blue if padj<0.01, red if log2FC>1 and padj<0.05)
with(subset(res, padj<.05 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(res, padj<.05 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))

#vsdata <- vst(dds, blind=FALSE)
#plotPCA(vsdata, intgroup="Cervix.Vagina") #using the DESEQ2 plotPCA fxn we can


#########################################################
# Compare the relative kegg expression
# between the cervical and vaginal samples using Deseq2
# per specificed CST
#########################################################


#read in kegg counts per sample
#mmtm_mt_comp <- read.csv("/Users/ichaudry/Library/CloudStorage/OneDrive-UniversityofMarylandSchoolofMedicine/research/ravel_lab/MMTM/analysis/01_mt_microbe_analysis/virgo2/taxa_filt_mmtm_mt_virgo2_kegg_counts_05132024.csv", header=TRUE, row.names=1)
mmtm_mt_comp <- read.csv("/Users/ichaudry/Library/CloudStorage/OneDrive-UniversityofMarylandSchoolofMedicine/research/ravel_lab/MMTM/analysis/01_mt_microbe_analysis/virgo2/taxa_filt_mmtm_mt_virgo2_kegg_counts_06122024_junk-genes-removed.csv", header=TRUE, row.names=1)


#filter samples with too few
min_reads <- 1000000
mmtm_mt_comp <- mmtm_mt_comp[,colSums(mmtm_mt_comp) >= min_reads]

#drop control columns
mmtm_mt_comp <- select(mmtm_mt_comp, starts_with("MT"))

#remove the MT_ and _ZYMO suffix from the column names
col_names <- colnames(mmtm_mt_comp)
col_names <- gsub("^MT_([A-Za-z0-9]+)_ZYMO$", "\\1", col_names)
colnames(mmtm_mt_comp) <- col_names

#load in the cst assingments per sample
#sample_csts <- read.csv('/Users/ichaudry/Documents/OneDrive - University of Maryland School of Medicine/research/ravel_lab/MMTM/analysis/cst_analysis/mmtm_cst_assignments_taxa_corr.csv',sep=',',row.names = 1, header = TRUE)
sample_csts <- read.csv('/Users/ichaudry/Documents/OneDrive - University of Maryland School of Medicine/research/ravel_lab/MMTM/analysis/01_mt_microbe_analysis/cst_assignments/mmtm_mt_cst_assingments_06122024_renamed.csv',sep=',',row.names = 1, header = TRUE)


#load the metadata for each sample
mmtm_metadata <- read.csv("/Users/ichaudry/Documents/OneDrive - University of Maryland School of Medicine/research/ravel_lab/MMTM/analysis/mmtm_metadata.csv", sep=',', row.names = 1)

#filter out the metadata dataframe to only include the samples that have composition data
samples = colnames(mmtm_mt_comp)
mmtm_metadata <- mmtm_metadata[row.names(mmtm_metadata) %in% samples, , drop=FALSE]
mmtm_metadata$Cervix.Vagina <- as.factor(mmtm_metadata$Cervix.Vagina)
mmtm_metadata$Timepoint <- as.factor(mmtm_metadata$Timepoint)
mmtm_metadata$Ext.Participant.ID <- as.factor(mmtm_metadata$Ext.Participant.ID)

#merge the CST information with the metadata
mmtm_metadata <- merge(sample_csts, mmtm_metadata, by='row.names' )
rownames(mmtm_metadata) <- mmtm_metadata$Row.names
mmtm_metadata$CST <- as.factor(mmtm_metadata$CST)

#filter out the control samples from the composition data by only keeping the samples that also have metadata associated with them
samples = rownames(mmtm_metadata)
mmtm_mt_comp <- mmtm_mt_comp[, colnames(mmtm_mt_comp) %in% samples, drop=FALSE]

#drop keggs with low abundance
mmtm_mt_comp <- mmtm_mt_comp[rowSums(mmtm_mt_comp[]) > 50,]


target_cst = c('IV-A', 'IV-B', 'IV-C')
mmtm_metadata <- mmtm_metadata[mmtm_metadata$CST %in% target_cst,]
samples = rownames(mmtm_metadata)
mmtm_mt_comp <- mmtm_mt_comp[, colnames(mmtm_mt_comp) %in% samples, drop=FALSE]


mmtm_mt_comp <- as.matrix(mmtm_mt_comp)+1
#mmtm_metadata <- as.matrix(mmtm_metadata)

dds <- DESeqDataSetFromMatrix(countData=mmtm_mt_comp, 
                              colData=mmtm_metadata, 
                              design=~Timepoint + Cervix.Vagina)

dds <- DESeq(dds)



res <- results(dds, contrast = c("Cervix.Vagina", "C", "V"), tidy = TRUE)
head(res) #let's look at the results table


res <- res[order(res$padj),]
head(res)

write.csv(res, '/Users/ichaudry/Library/CloudStorage/OneDrive-UniversityofMarylandSchoolofMedicine/research/ravel_lab/MMTM/analysis/01_mt_microbe_analysis/deseq2/mmtm_mt_diff_micro_kegg_cervix-vagina_cstIV_deseq2.csv', row.names = FALSE)

#reset par
par(mfrow=c(1,1))
# Make a basic volcano plot
with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-5,5)))

# Add colored points: blue if padj<0.01, red if log2FC>1 and padj<0.05)
with(subset(res, padj<.05 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(res, padj<.05 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))

##



#############################################################
# Identify differentially expressed genes per taxa 
# between the cervix and vagina. 
#############################################################

completed_taxa <- c()
#load the metagenomics gene level abundance data determined by mapping reads to VIRGO2
#mmtm_mt_comp_all <- read.csv("/Users/ichaudry/Library/CloudStorage/OneDrive-UniversityofMarylandSchoolofMedicine/research/ravel_lab/MMTM/analysis/01_mt_microbe_analysis/virgo2/taxa_filt_mmtm_mt_virgo2_NR_Taxa_051324.csv", header=TRUE, row.names=1)
mmtm_mt_comp_all <- read.csv("/Users/ichaudry/Library/CloudStorage/OneDrive-UniversityofMarylandSchoolofMedicine/research/ravel_lab/MMTM/analysis/01_mt_microbe_analysis/virgo2/taxa_filt_mmtm_mt_virgo2_NR_Taxa_051324_junk_removed.csv", header=TRUE, row.names=1)


#drop all the genes with unassigned taxa 
mmtm_mt_comp_all <- subset(mmtm_mt_comp_all, Taxa != "")

#top taxa to be used for per-taxa DGE analysi
target_taxa <- c('Lactobacillus_iners', 'MultiGenera', 'Megasphaera_lornae', 'UBA629_sp005465875', 'Gardnerella_vaginalis', 'Prevotella_amnii', 'Prevotella_spNov3', 'Gardnerella_swidsinkii', 'Lactobacillus_crispatus', 'Gardnerella_vaginalis_A', 'Roseburia_hominis;Roseburia', 'Gardnerella_vaginalis_D', 'Gardnerella', 'Gardnerella_leopoldii', 'Gardnerella_vaginalis_C', 'Gardnerella_vaginalis_H', 'Gardnerella_vaginalis_F', 'Prevotella_timonensis')
#target_taxa <- c('Lactobacillus_crispatus', 'Gardnerella_leopoldii')

#setting the minimum number of reads for the specified taxa needed for the sample to be considered
min_reads = 500000

#iterate through all the taxa and do the DGE for each one
for(taxa in target_taxa){ 
  print(paste("Starting ", taxa, sep=""))
  
  if (taxa %in% completed_taxa){
    next
  }
  
  #filter out genes only from this taxa
  mmtm_mt_comp <- subset(mmtm_mt_comp_all, Taxa == taxa)
  #mmtm_mt_comp <- subset(mmtm_mt_comp_all, grepl(paste0("^", taxa), Taxa))
  
  #drop the category and length columns from the mapping table
  mmtm_mt_comp <- subset(mmtm_mt_comp, select=-c(Taxa, Cat, Length))
  
  #drop genes that are mostly 0s
  #mmtm_mt_comp <- mmtm_mt_comp[(rowSums(mmtm_mt_comp > 0) / ncol(mmtm_mt_comp) > 0.15), ]
  
  #drop samples that do not meet the min_read requirement
  #mmtm_mt_comp <- mmtm_mt_comp[, !apply(mmtm_mt_comp == 0, 2, all)]
  mmtm_mt_comp <- mmtm_mt_comp[,colSums(mmtm_mt_comp) >= min_reads]
  
  #drop control columns
  mmtm_mt_comp <- select(mmtm_mt_comp, starts_with("MT"))
  
  #remove the MT_ and _ZYMO suffix from the column names
  col_names <- colnames(mmtm_mt_comp)
  col_names <- gsub("^MT_([A-Za-z0-9]+)_ZYMO$", "\\1", col_names)
  colnames(mmtm_mt_comp) <- col_names
  
  #load in the cst assingments per sample
  #sample_csts <- read.csv('/Users/ichaudry/Documents/OneDrive - University of Maryland School of Medicine/research/ravel_lab/MMTM/analysis/cst_analysis/mmtm_cst_assignments_taxa_corr.csv',sep=',',row.names = 1, header = TRUE)
  sample_csts <- read.csv('/Users/ichaudry/Documents/OneDrive - University of Maryland School of Medicine/research/ravel_lab/MMTM/analysis/01_mt_microbe_analysis/cst_assignments/mmtm_mt_cst_assingments_06122024_renamed.csv',sep=',',row.names = 1, header = TRUE)
  
  #load the metadata for each sample
  mmtm_metadata <- read.csv("/Users/ichaudry/Documents/OneDrive - University of Maryland School of Medicine/research/ravel_lab/MMTM/analysis/mmtm_metadata.csv", sep=',', row.names = 1)
  
  #filter out the metadata dataframe to only include the samples that have composition data
  samples = colnames(mmtm_mt_comp)
  mmtm_metadata <- mmtm_metadata[row.names(mmtm_metadata) %in% samples, , drop=FALSE]
  mmtm_metadata$Cervix.Vagina <- as.factor(mmtm_metadata$Cervix.Vagina)
  mmtm_metadata$Timepoint <- as.factor(mmtm_metadata$Timepoint)
  mmtm_metadata$Ext.Participant.ID <- as.factor(mmtm_metadata$Ext.Participant.ID)
  
  #merge the CST information with the metadata
  mmtm_metadata <- merge(sample_csts, mmtm_metadata, by='row.names' )
  rownames(mmtm_metadata) <- mmtm_metadata$Row.names
  mmtm_metadata$CST <- as.factor(mmtm_metadata$CST)
  
  #filter out the control samples from the composition data by only keeping the samples that also have metadata associated with them
  samples = rownames(mmtm_metadata)
  mmtm_mt_comp <- mmtm_mt_comp[, colnames(mmtm_mt_comp) %in% samples, drop=FALSE]
  
  mmtm_mt_comp <- as.matrix(mmtm_mt_comp)+1
  #dds <- DESeqDataSetFromMatrix(countData=mmtm_mt_comp, 
  #                              colData=mmtm_metadata, 
  #                              design=~CST + Timepoint + Cervix.Vagina)
  
  dds <- DESeqDataSetFromMatrix(countData=mmtm_mt_comp, 
                                colData=mmtm_metadata, 
                                design=~Timepoint + Cervix.Vagina)
    
  dds <- DESeq(dds)
    
    
    
  res <- results(dds, contrast = c("Cervix.Vagina", "C", "V"), tidy = TRUE)
  head(res) #let's look at the results table
  
  
  res <- res[order(res$padj),]
  head(res)
  
  write.csv(res, paste('/Users/ichaudry/Library/CloudStorage/OneDrive-UniversityofMarylandSchoolofMedicine/research/ravel_lab/MMTM/analysis/01_mt_microbe_analysis/deseq2/per_taxa/mmtm_mt_diff_micro_genes_cervix-vagina_deseq2_', taxa, '.csv', sep = ""), row.names = FALSE)
  
  #reset par
  par(mfrow=c(1,1))
  
  # Make a basic volcano plot
  #with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main=taxa, xlim=c(-5,5)))
  
  # Add colored points: blue if padj<0.01, red if log2FC>1 and padj<0.05)
  #with(subset(res, padj<.05 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
  #with(subset(res, padj<.05 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
  
  #vsdata <- vst(dds, blind=FALSE)
  #plotPCA(vsdata, intgroup="Cervix.Vagina") #using the DESEQ2 plotPCA fxn we can
  
  completed_taxa <- c(completed_taxa, taxa)
}


#########################################################
#########################################################
# OVERTIME ANALYSIS
#########################################################
#########################################################

#########################################################
# Compare the relative kegg expression
# in the cervical and vaginal samples OVERTIME using Deseq2
# looking across all CSTs, but a single location
#########################################################

#read in kegg counts per sample
#mmtm_mt_comp <- read.csv("/Users/ichaudry/Library/CloudStorage/OneDrive-UniversityofMarylandSchoolofMedicine/research/ravel_lab/MMTM/analysis/01_mt_microbe_analysis/virgo2/taxa_filt_mmtm_mt_virgo2_kegg_counts_05132024.csv", header=TRUE, row.names=1)
mmtm_mt_comp <- read.csv("/Users/ichaudry/Library/CloudStorage/OneDrive-UniversityofMarylandSchoolofMedicine/research/ravel_lab/MMTM/analysis/01_mt_microbe_analysis/virgo2/taxa_filt_mmtm_mt_virgo2_kegg_counts_06122024_junk-genes-removed.csv", header=TRUE, row.names=1)


#filter samples with too few
min_reads <- 1000000
mmtm_mt_comp <- mmtm_mt_comp[,colSums(mmtm_mt_comp) >= min_reads]

#drop control columns
mmtm_mt_comp <- select(mmtm_mt_comp, starts_with("MT"))

#remove the MT_ and _ZYMO suffix from the column names
col_names <- colnames(mmtm_mt_comp)
col_names <- gsub("^MT_([A-Za-z0-9]+)_ZYMO$", "\\1", col_names)
colnames(mmtm_mt_comp) <- col_names

#load in the cst assingments per sample
#sample_csts <- read.csv('/Users/ichaudry/Documents/OneDrive - University of Maryland School of Medicine/research/ravel_lab/MMTM/analysis/cst_analysis/mmtm_cst_assignments_taxa_corr.csv',sep=',',row.names = 1, header = TRUE)
sample_csts <- read.csv('/Users/ichaudry/Documents/OneDrive - University of Maryland School of Medicine/research/ravel_lab/MMTM/analysis/01_mt_microbe_analysis/cst_assignments/mmtm_mt_cst_assingments_06122024_renamed.csv',sep=',',row.names = 1, header = TRUE)


#load the metadata for each sample
mmtm_metadata <- read.csv("/Users/ichaudry/Documents/OneDrive - University of Maryland School of Medicine/research/ravel_lab/MMTM/analysis/mmtm_metadata.csv", sep=',', row.names = 1)

#filter out the metadata dataframe to only include the samples that have composition data
samples = colnames(mmtm_mt_comp)
mmtm_metadata <- mmtm_metadata[row.names(mmtm_metadata) %in% samples, , drop=FALSE]
mmtm_metadata$Cervix.Vagina <- as.factor(mmtm_metadata$Cervix.Vagina)
mmtm_metadata$Timepoint <- as.factor(mmtm_metadata$Timepoint)
mmtm_metadata$Ext.Participant.ID <- as.factor(mmtm_metadata$Ext.Participant.ID)

#merge the CST information with the metadata
mmtm_metadata <- merge(sample_csts, mmtm_metadata, by='row.names' )
rownames(mmtm_metadata) <- mmtm_metadata$Row.names
mmtm_metadata$CST <- as.factor(mmtm_metadata$CST)

#filter out the control samples from the composition data by only keeping the samples that also have metadata associated with them
samples = rownames(mmtm_metadata)
mmtm_mt_comp <- mmtm_mt_comp[, colnames(mmtm_mt_comp) %in% samples, drop=FALSE]

#drop keggs with low abundance
mmtm_mt_comp <- mmtm_mt_comp[rowSums(mmtm_mt_comp[]) > 50,]


target_location = c('V')
mmtm_metadata <- mmtm_metadata[mmtm_metadata$Cervix.Vagina %in% target_location,]
samples = rownames(mmtm_metadata)
mmtm_mt_comp <- mmtm_mt_comp[, colnames(mmtm_mt_comp) %in% samples, drop=FALSE]


mmtm_mt_comp <- as.matrix(mmtm_mt_comp)+1
#mmtm_metadata <- as.matrix(mmtm_metadata)

dds <- DESeqDataSetFromMatrix(countData=mmtm_mt_comp, 
                              colData=mmtm_metadata, 
                              design=~CST+Timepoint)

dds <- DESeq(dds)



res <- results(dds, contrast = c("Timepoint", "1", "2"), tidy = TRUE)
head(res) #let's look at the results table


res <- res[order(res$padj),]
head(res)

write.csv(res, '/Users/ichaudry/Library/CloudStorage/OneDrive-UniversityofMarylandSchoolofMedicine/research/ravel_lab/MMTM/analysis/01_mt_microbe_analysis/deseq2/mmtm_mt_diff_micro_kegg_overtime_vagina_deseq2.csv', row.names = FALSE)

#reset par
par(mfrow=c(1,1))
# Make a basic volcano plot
with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-5,5)))

# Add colored points: blue if padj<0.01, red if log2FC>1 and padj<0.05)
with(subset(res, padj<.05 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(res, padj<.05 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))

##


#########################################################
# Compare the relative kegg expression
# in the cervical and vaginal samples OVERTIME using Deseq2
# per specified CST
#########################################################

#read in kegg counts per sample
#mmtm_mt_comp <- read.csv("/Users/ichaudry/Library/CloudStorage/OneDrive-UniversityofMarylandSchoolofMedicine/research/ravel_lab/MMTM/analysis/01_mt_microbe_analysis/virgo2/taxa_filt_mmtm_mt_virgo2_kegg_counts_05132024.csv", header=TRUE, row.names=1)
mmtm_mt_comp <- read.csv("/Users/ichaudry/Library/CloudStorage/OneDrive-UniversityofMarylandSchoolofMedicine/research/ravel_lab/MMTM/analysis/01_mt_microbe_analysis/virgo2/taxa_filt_mmtm_mt_virgo2_kegg_counts_06122024_junk-genes-removed.csv", header=TRUE, row.names=1)


#filter samples with too few
min_reads <- 1000000
mmtm_mt_comp <- mmtm_mt_comp[,colSums(mmtm_mt_comp) >= min_reads]

#drop control columns
mmtm_mt_comp <- select(mmtm_mt_comp, starts_with("MT"))

#remove the MT_ and _ZYMO suffix from the column names
col_names <- colnames(mmtm_mt_comp)
col_names <- gsub("^MT_([A-Za-z0-9]+)_ZYMO$", "\\1", col_names)
colnames(mmtm_mt_comp) <- col_names

#load in the cst assingments per sample
#sample_csts <- read.csv('/Users/ichaudry/Documents/OneDrive - University of Maryland School of Medicine/research/ravel_lab/MMTM/analysis/cst_analysis/mmtm_cst_assignments_taxa_corr.csv',sep=',',row.names = 1, header = TRUE)
sample_csts <- read.csv('/Users/ichaudry/Documents/OneDrive - University of Maryland School of Medicine/research/ravel_lab/MMTM/analysis/01_mt_microbe_analysis/cst_assignments/mmtm_mt_cst_assingments_06122024_renamed.csv',sep=',',row.names = 1, header = TRUE)


#load the metadata for each sample
mmtm_metadata <- read.csv("/Users/ichaudry/Documents/OneDrive - University of Maryland School of Medicine/research/ravel_lab/MMTM/analysis/mmtm_metadata.csv", sep=',', row.names = 1)

#filter out the metadata dataframe to only include the samples that have composition data
samples = colnames(mmtm_mt_comp)
mmtm_metadata <- mmtm_metadata[row.names(mmtm_metadata) %in% samples, , drop=FALSE]
mmtm_metadata$Cervix.Vagina <- as.factor(mmtm_metadata$Cervix.Vagina)
mmtm_metadata$Timepoint <- as.factor(mmtm_metadata$Timepoint)
mmtm_metadata$Ext.Participant.ID <- as.factor(mmtm_metadata$Ext.Participant.ID)

#merge the CST information with the metadata
mmtm_metadata <- merge(sample_csts, mmtm_metadata, by='row.names' )
rownames(mmtm_metadata) <- mmtm_metadata$Row.names
mmtm_metadata$CST <- as.factor(mmtm_metadata$CST)

#filter out the control samples from the composition data by only keeping the samples that also have metadata associated with them
samples = rownames(mmtm_metadata)
mmtm_mt_comp <- mmtm_mt_comp[, colnames(mmtm_mt_comp) %in% samples, drop=FALSE]

#drop keggs with low abundance
mmtm_mt_comp <- mmtm_mt_comp[rowSums(mmtm_mt_comp[]) > 50,]


target_cst = c('IV-A', 'IV-B', 'IV-C')
mmtm_metadata <- mmtm_metadata[mmtm_metadata$CST %in% target_cst,]
samples = rownames(mmtm_metadata)
mmtm_mt_comp <- mmtm_mt_comp[, colnames(mmtm_mt_comp) %in% samples, drop=FALSE]

mmtm_mt_comp <- as.matrix(mmtm_mt_comp)+1
#mmtm_metadata <- as.matrix(mmtm_metadata)

dds <- DESeqDataSetFromMatrix(countData=mmtm_mt_comp, 
                              colData=mmtm_metadata, 
                              design=~Cervix.Vagina+Timepoint)

dds <- DESeq(dds)



res <- results(dds, contrast = c("Timepoint", "1", "2"), tidy = TRUE)
head(res) #let's look at the results table


res <- res[order(res$padj),]
head(res)

write.csv(res, '/Users/ichaudry/Library/CloudStorage/OneDrive-UniversityofMarylandSchoolofMedicine/research/ravel_lab/MMTM/analysis/01_mt_microbe_analysis/deseq2/mmtm_mt_diff_micro_kegg_overtime_cstIV_deseq2.csv', row.names = FALSE)

#reset par
par(mfrow=c(1,1))
# Make a basic volcano plot
with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-5,5)))

# Add colored points: blue if padj<0.01, red if log2FC>1 and padj<0.05)
with(subset(res, padj<.05 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(res, padj<.05 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))


#########################################################
# Identify differentially expressed genes per taxa 
# between the timepoints 1  and 2. 
#########################################################

completed_taxa <- c()
#load the metagenomics gene level abundance data determined by mapping reads to VIRGO2
#mmtm_mt_comp_all <- read.csv("/Users/ichaudry/Library/CloudStorage/OneDrive-UniversityofMarylandSchoolofMedicine/research/ravel_lab/MMTM/analysis/01_mt_microbe_analysis/virgo2/taxa_filt_mmtm_mt_virgo2_NR_Taxa_051324.csv", header=TRUE, row.names=1)
mmtm_mt_comp_all <- read.csv("/Users/ichaudry/Library/CloudStorage/OneDrive-UniversityofMarylandSchoolofMedicine/research/ravel_lab/MMTM/analysis/01_mt_microbe_analysis/virgo2/taxa_filt_mmtm_mt_virgo2_NR_Taxa_051324_junk_removed.csv", header=TRUE, row.names=1)


#drop all the genes with unassigned taxa 
mmtm_mt_comp_all <- subset(mmtm_mt_comp_all, Taxa != "")

#top taxa to be used for per-taxa DGE analysi
target_taxa <- c('Lactobacillus_iners', 'MultiGenera', 'Megasphaera_lornae', 'UBA629_sp005465875', 'Gardnerella_vaginalis', 'Prevotella_amnii', 'Prevotella_spNov3', 'Gardnerella_swidsinkii', 'Lactobacillus_crispatus', 'Gardnerella_vaginalis_A', 'Roseburia_hominis;Roseburia', 'Gardnerella_vaginalis_D', 'Gardnerella', 'Gardnerella_leopoldii', 'Gardnerella_vaginalis_C', 'Gardnerella_vaginalis_H', 'Gardnerella_vaginalis_F', 'Prevotella_timonensis')
#target_taxa <- c('Lactobacillus_crispatus', 'Gardnerella_leopoldii')

#setting the minimum number of reads for the specified taxa needed for the sample to be considered
min_reads = 500000

#iterate through all the taxa and do the DGE for each one
for(taxa in target_taxa){ 
  print(paste("Starting ", taxa, sep=""))
  
  if (taxa %in% completed_taxa){
    next
  }
  
  #filter out genes only from this taxa
  mmtm_mt_comp <- subset(mmtm_mt_comp_all, Taxa == taxa)
  #mmtm_mt_comp <- subset(mmtm_mt_comp_all, grepl(paste0("^", taxa), Taxa))
  
  #drop the category and length columns from the mapping table
  mmtm_mt_comp <- subset(mmtm_mt_comp, select=-c(Taxa, Cat, Length))
  
  #drop genes that are mostly 0s
  #mmtm_mt_comp <- mmtm_mt_comp[(rowSums(mmtm_mt_comp > 0) / ncol(mmtm_mt_comp) > 0.15), ]
  
  #drop samples that do not meet the min_read requirement
  #mmtm_mt_comp <- mmtm_mt_comp[, !apply(mmtm_mt_comp == 0, 2, all)]
  mmtm_mt_comp <- mmtm_mt_comp[,colSums(mmtm_mt_comp) >= min_reads]
  
  #drop control columns
  mmtm_mt_comp <- select(mmtm_mt_comp, starts_with("MT"))
  
  #remove the MT_ and _ZYMO suffix from the column names
  col_names <- colnames(mmtm_mt_comp)
  col_names <- gsub("^MT_([A-Za-z0-9]+)_ZYMO$", "\\1", col_names)
  colnames(mmtm_mt_comp) <- col_names
  
  #load in the cst assingments per sample
  #sample_csts <- read.csv('/Users/ichaudry/Documents/OneDrive - University of Maryland School of Medicine/research/ravel_lab/MMTM/analysis/cst_analysis/mmtm_cst_assignments_taxa_corr.csv',sep=',',row.names = 1, header = TRUE)
  sample_csts <- read.csv('/Users/ichaudry/Documents/OneDrive - University of Maryland School of Medicine/research/ravel_lab/MMTM/analysis/01_mt_microbe_analysis/cst_assignments/mmtm_mt_cst_assingments_06122024_renamed.csv',sep=',',row.names = 1, header = TRUE)
  
  #load the metadata for each sample
  mmtm_metadata <- read.csv("/Users/ichaudry/Documents/OneDrive - University of Maryland School of Medicine/research/ravel_lab/MMTM/analysis/mmtm_metadata.csv", sep=',', row.names = 1)
  
  #filter out the metadata dataframe to only include the samples that have composition data
  samples = colnames(mmtm_mt_comp)
  mmtm_metadata <- mmtm_metadata[row.names(mmtm_metadata) %in% samples, , drop=FALSE]
  mmtm_metadata$Cervix.Vagina <- as.factor(mmtm_metadata$Cervix.Vagina)
  mmtm_metadata$Timepoint <- as.factor(mmtm_metadata$Timepoint)
  mmtm_metadata$Ext.Participant.ID <- as.factor(mmtm_metadata$Ext.Participant.ID)
  
  #merge the CST information with the metadata
  mmtm_metadata <- merge(sample_csts, mmtm_metadata, by='row.names' )
  rownames(mmtm_metadata) <- mmtm_metadata$Row.names
  mmtm_metadata$CST <- as.factor(mmtm_metadata$CST)
  
  #filter out the control samples from the composition data by only keeping the samples that also have metadata associated with them
  samples = rownames(mmtm_metadata)
  mmtm_mt_comp <- mmtm_mt_comp[, colnames(mmtm_mt_comp) %in% samples, drop=FALSE]
  
  mmtm_mt_comp <- as.matrix(mmtm_mt_comp)+1
  #dds <- DESeqDataSetFromMatrix(countData=mmtm_mt_comp, 
  #                              colData=mmtm_metadata, 
  #                              design=~CST + Timepoint + Cervix.Vagina)
  
  dds <- DESeqDataSetFromMatrix(countData=mmtm_mt_comp, 
                                colData=mmtm_metadata, 
                                design=~Cervix.Vagina + Timepoint)
  
  dds <- DESeq(dds)
  
  
  
  res <- results(dds, contrast = c("Timepoint", "1", "2"), tidy = TRUE)
  head(res) #let's look at the results table
  
  
  res <- res[order(res$padj),]
  head(res)
  
  write.csv(res, paste('/Users/ichaudry/Library/CloudStorage/OneDrive-UniversityofMarylandSchoolofMedicine/research/ravel_lab/MMTM/analysis/01_mt_microbe_analysis/deseq2/per_taxa/mmtm_mt_diff_micro_genes_overtime_deseq2_', taxa, '.csv', sep = ""), row.names = FALSE)
  
  #reset par
  par(mfrow=c(1,1))
  
  # Make a basic volcano plot
  #with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main=taxa, xlim=c(-5,5)))
  
  # Add colored points: blue if padj<0.01, red if log2FC>1 and padj<0.05)
  #with(subset(res, padj<.05 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
  #with(subset(res, padj<.05 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
  
  #vsdata <- vst(dds, blind=FALSE)
  #plotPCA(vsdata, intgroup="Cervix.Vagina") #using the DESEQ2 plotPCA fxn we can
  
  completed_taxa <- c(completed_taxa, taxa)
}



