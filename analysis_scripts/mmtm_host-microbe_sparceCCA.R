library("PMA")
library("kernlab")
library(dplyr)
library(compositions)
setwd('/Users/ichaudry/Documents/OneDrive - University of Maryland School of Medicine/research/ravel_lab/MMTM/analysis/04_host_microbe_mt_integ')


# read in sample metadata
########################
sample_csts <- read.csv('/Users/ichaudry/Documents/OneDrive - University of Maryland School of Medicine/research/ravel_lab/MMTM/analysis/01_mt_microbe_analysis/cst_assignments/mmtm_mt_cst_assingments_06122024_renamed.csv',sep=',', header = TRUE)
row.names(sample_csts) <- sample_csts$sampleID

#load the metadata for each sample
mmtm_metadata <- read.csv("/Users/ichaudry/Documents/OneDrive - University of Maryland School of Medicine/research/ravel_lab/MMTM/analysis/mmtm_metadata.csv")
rownames(mmtm_metadata) <- mmtm_metadata$IGS.LABEL
mmtm_metadata <- mmtm_metadata[,-1]


# host transcriptomics and microbial relative abundance
#######################################################

# read in the human transcriptomics
mmtm_mt_host <- read.csv("/Users/ichaudry/Documents/OneDrive - University of Maryland School of Medicine/research/ravel_lab/MMTM/analysis/02_mt_host_analysis/MMTM_host_featureCounts_0520204_tpms.csv")
rownames(mmtm_mt_host) <- mmtm_mt_host$Geneid
mmtm_mt_host <- mmtm_mt_host[,-1]
mmtm_mt_host <- mmtm_mt_host[rowSums(mmtm_mt_host[]) > 0,]
mmtm_mt_host <- mmtm_mt_host[rowMeans(mmtm_mt_host[]) > 0.5,]
mmtm_mt_host <- mmtm_mt_host %>% select(starts_with('MTM'))
mmtm_mt_host <- data.frame(t(mmtm_mt_host))

#read in the microbial relative abundance 
mmtm_mg_comp <- read.csv("/Users/ichaudry/Documents/OneDrive - University of Maryland School of Medicine/research/ravel_lab/MMTM/analysis/00_mg_analysis/virgo2/taxa_filt_mmtm_mg_virgo2.csv")
rownames(mmtm_mg_comp) <- mmtm_mg_comp$Taxa
mmtm_mg_comp <- mmtm_mg_comp[,-1]
mmtm_mg_comp <- mmtm_mg_comp[rowSums(mmtm_mg_comp[]) > 0,]
mmtm_mg_comp <- mmtm_mg_comp[rowMeans(mmtm_mg_comp[]) > 1E-3,]
mmtm_mg_comp <- mmtm_mg_comp %>% select(starts_with('MTM'))
mmtm_mg_comp <- data.frame(t(mmtm_mg_comp))
mmtm_mg_comp <- clr(mmtm_mg_comp)
mmtm_mg_comp <- data.frame(mmtm_mg_comp)

# Find the intersection of the samples
common_rows <- intersect(rownames(mmtm_mt_host), rownames(mmtm_mg_comp))

# Subset to keep only the common rows
mmtm_mt_host <- mmtm_mt_host[common_rows, , drop = FALSE]
mmtm_mg_comp <- mmtm_mg_comp[common_rows, , drop = FALSE]

#filter by bodysite
target_samples <- mmtm_metadata[mmtm_metadata$Cervix.Vagina == 'V',]
target_samples <- rownames(target_samples)

mmtm_mt_host <- mmtm_mt_host[rownames(mmtm_mt_host) %in% target_samples, , drop = FALSE]
mmtm_mg_comp <- mmtm_mg_comp[rownames(mmtm_mg_comp) %in% target_samples, , drop = FALSE]

#set up dataframes
x_data <- mmtm_mt_host
x_data[is.na(x_data)] <- 0

z_data <- mmtm_mg_comp
z_data[is.na(z_data)] <- 0

x_data <- x_data[, apply(x_data, 2, var) != 0]
z_data <- z_data[, apply(z_data, 2, var) != 0]


#check for association between the two datasets
cca_permutation <- CCA.permute(x=as.matrix(x_data),z=as.matrix(z_data),typex = "standard",typez = "standard",niter=3,trace=TRUE,standardize=TRUE,
                               penaltyxs =c(0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.075,0.075,0.075,0.075,0.075,0.075,0.075,0.1,0.1,0.1,0.1,0.1,0.1,.1,0.125,0.125,0.125,0.125,0.125,0.125,.125,0.15,0.15,0.15,0.15,0.15,0.15,.15,.2,.2,.2,.2,.2,.2,.2),
                               penaltyzs = c(0.1, 0.2,0.3,0.4,0.5,0.6,0.7, 0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.1,0.2,0.3,0.4,0.5,0.6,0.7))
plot(cca_permutation)
cca_permutation

cca_output <- CCA(x=x_data,z=z_data,typex = "standard",typez = "standard",niter=3,trace=TRUE,standardize=TRUE, K=5, 
                  penaltyx=cca_permutation$bestpenaltyx, penaltyz=cca_permutation$bestpenaltyz)

cca_output <- CCA(x=x_data,z=z_data,typex = "standard",typez = "standard",niter=3,trace=TRUE,standardize=TRUE, K=5, 
                  penaltyx=0.1, penaltyz=0.5)


u_vals <- data.frame(cca_output$u)
rownames(u_vals) <- cca_output$xnames

v_vals <- data.frame(cca_output$v)
rownames(v_vals) <- cca_output$znames







# host transcriptomics and microbial kegg expression
#######################################################

# read in the human transcriptomics
mmtm_mt_host <- read.csv("/Users/ichaudry/Documents/OneDrive - University of Maryland School of Medicine/research/ravel_lab/MMTM/analysis/02_mt_host_analysis/MMTM_host_featureCounts_0520204_tpms.csv")
rownames(mmtm_mt_host) <- mmtm_mt_host$Geneid
mmtm_mt_host <- mmtm_mt_host[,-1]
mmtm_mt_host <- mmtm_mt_host[rowSums(mmtm_mt_host[]) > 0,]
mmtm_mt_host <- mmtm_mt_host[rowMeans(mmtm_mt_host[]) > 1,]
mmtm_mt_host <- mmtm_mt_host %>% select(starts_with('MTM'))
mmtm_mt_host <- data.frame(t(mmtm_mt_host))

#read in the microbial kegg expression 
mmtm_mt_kegg_tpms <- read.csv("/Users/ichaudry/Library/CloudStorage/OneDrive-UniversityofMarylandSchoolofMedicine/research/ravel_lab/MMTM/analysis/01_mt_microbe_analysis/virgo2/taxa_filt_mmtm_mt_virgo2_kegg_tpms_06122024_junk-genes-removed.csv", header=TRUE, row.names=1)
mmtm_mt_kegg_tpms <- mmtm_mt_kegg_tpms[rowMeans(mmtm_mt_kegg_tpms[]) > 1,]
mmtm_mt_kegg_tpms <- mmtm_mt_kegg_tpms %>% select(starts_with('MT_MTM'))
colnames(mmtm_mt_kegg_tpms) <- gsub("MT_|_ZYMO", "", colnames(mmtm_mt_kegg_tpms))
mmtm_mt_kegg_tpms <- data.frame(t(mmtm_mt_kegg_tpms))

# Find the intersection of the samples
common_rows <- intersect(rownames(mmtm_mt_host), rownames(mmtm_mt_kegg_tpms))

# Subset to keep only the common rows
mmtm_mt_host <- mmtm_mt_host[common_rows, , drop = FALSE]
mmtm_mt_kegg_tpms <- mmtm_mt_kegg_tpms[common_rows, , drop = FALSE]

#set up dataframes
x_data <- mmtm_mt_host
x_data[is.na(x_data)] <- 0

z_data <- mmtm_mt_kegg_tpms
z_data[is.na(z_data)] <- 0

x_data <- x_data[, apply(x_data, 2, var) != 0]
z_data <- z_data[, apply(z_data, 2, var) != 0]



#check for association between the two datasets
cca_permutation <- CCA.permute(x=as.matrix(x_data),z=as.matrix(z_data),typex = "standard",typez = "standard",niter=3,trace=TRUE,standardize=TRUE,
                               penaltyxs =c(0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.075,0.075,0.075,0.075,0.075,0.075,0.075,0.1,0.1,0.1,0.1,0.1,0.1,.1,0.125,0.125,0.125,0.125,0.125,0.125,.125,0.15,0.15,0.15,0.15,0.15,0.15,.15,.2,.2,.2,.2,.2,.2,.2),
                               penaltyzs = c(0.1, 0.2,0.3,0.4,0.5,0.6,0.7, 0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.1,0.2,0.3,0.4,0.5,0.6,0.7))
plot(cca_permutation)
cca_permutation

cca_output <- CCA(x=x_data,z=z_data,typex = "standard",typez = "standard",niter=3,trace=TRUE,standardize=TRUE, K=5, 
                  penaltyx=cca_permutation$bestpenaltyx, penaltyz=cca_permutation$bestpenaltyz)


u_vals <- data.frame(cca_output$u)
rownames(u_vals) <- cca_output$xnames

v_vals <- data.frame(cca_output$v)
rownames(v_vals) <- cca_output$znames

cca_output

write.csv(u_vals, "host_mt_gene-microbe_mt_kegg_HOST.csv")
write.csv(v_vals, "host_mt_gene-microbe_mt_kegg_MICROBE.csv")
write.csv(cca_output$cors,'host_mt_gene-microbe_mt_kegg_CORR.csv')




# host transcriptomics and microbial kegg expression: BODYSITE
#######################################################

# read in the human transcriptomics
mmtm_mt_host <- read.csv("/Users/ichaudry/Documents/OneDrive - University of Maryland School of Medicine/research/ravel_lab/MMTM/analysis/02_mt_host_analysis/MMTM_host_featureCounts_0520204_tpms.csv")
rownames(mmtm_mt_host) <- mmtm_mt_host$Geneid
mmtm_mt_host <- mmtm_mt_host[,-1]
mmtm_mt_host <- mmtm_mt_host[rowSums(mmtm_mt_host[]) > 0,]
mmtm_mt_host <- mmtm_mt_host[rowMeans(mmtm_mt_host[]) > 1,]
mmtm_mt_host <- mmtm_mt_host %>% select(starts_with('MTM'))
mmtm_mt_host <- data.frame(t(mmtm_mt_host))

#read in the microbial kegg expression 
mmtm_mt_kegg_tpms <- read.csv("/Users/ichaudry/Library/CloudStorage/OneDrive-UniversityofMarylandSchoolofMedicine/research/ravel_lab/MMTM/analysis/01_mt_microbe_analysis/virgo2/taxa_filt_mmtm_mt_virgo2_kegg_tpms_06122024_junk-genes-removed.csv", header=TRUE, row.names=1)
mmtm_mt_kegg_tpms <- mmtm_mt_kegg_tpms[rowMeans(mmtm_mt_kegg_tpms[]) > 1,]
mmtm_mt_kegg_tpms <- mmtm_mt_kegg_tpms %>% select(starts_with('MT_MTM'))
colnames(mmtm_mt_kegg_tpms) <- gsub("MT_|_ZYMO", "", colnames(mmtm_mt_kegg_tpms))
mmtm_mt_kegg_tpms <- data.frame(t(mmtm_mt_kegg_tpms))

# Find the intersection of the samples
common_rows <- intersect(rownames(mmtm_mt_host), rownames(mmtm_mt_kegg_tpms))

# Subset to keep only the common rows
mmtm_mt_host <- mmtm_mt_host[common_rows, , drop = FALSE]
mmtm_mt_kegg_tpms <- mmtm_mt_kegg_tpms[common_rows, , drop = FALSE]

#filter by bodysite
target_samples <- mmtm_metadata[mmtm_metadata$Cervix.Vagina == 'C',]  ## <--- set body site here
target_samples <- rownames(target_samples)

mmtm_mt_host <- mmtm_mt_host[rownames(mmtm_mt_host) %in% target_samples, , drop = FALSE]
mmtm_mt_kegg_tpms <- mmtm_mt_kegg_tpms[rownames(mmtm_mt_kegg_tpms) %in% target_samples, , drop = FALSE]

#set up dataframes
x_data <- mmtm_mt_host
x_data[is.na(x_data)] <- 0

z_data <- mmtm_mt_kegg_tpms
z_data[is.na(z_data)] <- 0

x_data <- x_data[, apply(x_data, 2, var) != 0]
z_data <- z_data[, apply(z_data, 2, var) != 0]



#check for association between the two datasets
cca_permutation <- CCA.permute(x=as.matrix(x_data),z=as.matrix(z_data),typex = "standard",typez = "standard",niter=3,trace=TRUE,standardize=TRUE,
                               penaltyxs =c(0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.075,0.075,0.075,0.075,0.075,0.075,0.075,0.1,0.1,0.1,0.1,0.1,0.1,.1,0.125,0.125,0.125,0.125,0.125,0.125,.125,0.15,0.15,0.15,0.15,0.15,0.15,.15,.2,.2,.2,.2,.2,.2,.2),
                               penaltyzs = c(0.1, 0.2,0.3,0.4,0.5,0.6,0.7, 0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.1,0.2,0.3,0.4,0.5,0.6,0.7))
plot(cca_permutation)
cca_permutation

cca_output <- CCA(x=x_data,z=z_data,typex = "standard",typez = "standard",niter=3,trace=TRUE,standardize=TRUE, K=5, 
                  penaltyx=cca_permutation$bestpenaltyx, penaltyz=cca_permutation$bestpenaltyz)


u_vals <- data.frame(cca_output$u)
rownames(u_vals) <- cca_output$xnames

v_vals <- data.frame(cca_output$v)
rownames(v_vals) <- cca_output$znames

cca_output

write.csv(u_vals, "VAGINA_host_mt_gene-microbe_mt_kegg_HOST.csv")
write.csv(v_vals, "VAGINA_host_mt_gene-microbe_mt_kegg_MICROBE.csv")
write.csv(cca_output$cors,'VAGINA_host_mt_gene-microbe_mt_kegg_CORR.csv')



# host transcriptomics and microbial kegg expression: CST
#######################################################

# read in the human transcriptomics
mmtm_mt_host <- read.csv("/Users/ichaudry/Documents/OneDrive - University of Maryland School of Medicine/research/ravel_lab/MMTM/analysis/02_mt_host_analysis/MMTM_host_featureCounts_0520204_tpms.csv")
rownames(mmtm_mt_host) <- mmtm_mt_host$Geneid
mmtm_mt_host <- mmtm_mt_host[,-1]
mmtm_mt_host <- mmtm_mt_host[rowSums(mmtm_mt_host[]) > 0,]
mmtm_mt_host <- mmtm_mt_host[rowMeans(mmtm_mt_host[]) > 1,]
mmtm_mt_host <- mmtm_mt_host %>% select(starts_with('MTM'))
mmtm_mt_host <- data.frame(t(mmtm_mt_host))

#read in the microbial kegg expression 
mmtm_mt_kegg_tpms <- read.csv("/Users/ichaudry/Library/CloudStorage/OneDrive-UniversityofMarylandSchoolofMedicine/research/ravel_lab/MMTM/analysis/01_mt_microbe_analysis/virgo2/taxa_filt_mmtm_mt_virgo2_kegg_tpms_06122024_junk-genes-removed.csv", header=TRUE, row.names=1)
mmtm_mt_kegg_tpms <- mmtm_mt_kegg_tpms[rowMeans(mmtm_mt_kegg_tpms[]) > 1,]
mmtm_mt_kegg_tpms <- mmtm_mt_kegg_tpms %>% select(starts_with('MT_MTM'))
colnames(mmtm_mt_kegg_tpms) <- gsub("MT_|_ZYMO", "", colnames(mmtm_mt_kegg_tpms))
mmtm_mt_kegg_tpms <- data.frame(t(mmtm_mt_kegg_tpms))

# Find the intersection of the samples
common_rows <- intersect(rownames(mmtm_mt_host), rownames(mmtm_mt_kegg_tpms))

# Subset to keep only the common rows
mmtm_mt_host <- mmtm_mt_host[common_rows, , drop = FALSE]
mmtm_mt_kegg_tpms <- mmtm_mt_kegg_tpms[common_rows, , drop = FALSE]

#filter by CST
#csts = c('I', 'II', 'III', 'V') ## <--- set CSTs here
csts = c('IV-A', 'IV-B', 'IV-C') ## <--- set CSTs here
target_samples <- sample_csts[sample_csts$CST %in% csts, ]  
target_samples <- rownames(target_samples)

mmtm_mt_host <- mmtm_mt_host[rownames(mmtm_mt_host) %in% target_samples, , drop = FALSE]
mmtm_mt_kegg_tpms <- mmtm_mt_kegg_tpms[rownames(mmtm_mt_kegg_tpms) %in% target_samples, , drop = FALSE]

#set up dataframes
x_data <- mmtm_mt_host
x_data[is.na(x_data)] <- 0

z_data <- mmtm_mt_kegg_tpms
z_data[is.na(z_data)] <- 0

x_data <- x_data[, apply(x_data, 2, var) != 0]
z_data <- z_data[, apply(z_data, 2, var) != 0]



#check for association between the two datasets
cca_permutation <- CCA.permute(x=as.matrix(x_data),z=as.matrix(z_data),typex = "standard",typez = "standard",niter=3,trace=TRUE,standardize=TRUE,
                               penaltyxs =c(0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.075,0.075,0.075,0.075,0.075,0.075,0.075,0.1,0.1,0.1,0.1,0.1,0.1,.1,0.125,0.125,0.125,0.125,0.125,0.125,.125,0.15,0.15,0.15,0.15,0.15,0.15,.15,.2,.2,.2,.2,.2,.2,.2),
                               penaltyzs = c(0.1, 0.2,0.3,0.4,0.5,0.6,0.7, 0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.1,0.2,0.3,0.4,0.5,0.6,0.7))
plot(cca_permutation)
cca_permutation

cca_output <- CCA(x=x_data,z=z_data,typex = "standard",typez = "standard",niter=3,trace=TRUE,standardize=TRUE, K=5, 
                  penaltyx=cca_permutation$bestpenaltyx, penaltyz=cca_permutation$bestpenaltyz)


u_vals <- data.frame(cca_output$u)
rownames(u_vals) <- cca_output$xnames

v_vals <- data.frame(cca_output$v)
rownames(v_vals) <- cca_output$znames

cca_output

write.csv(u_vals, "CSTIV_host_mt_gene-microbe_mt_kegg_HOST.csv")
write.csv(v_vals, "CSTIV_host_mt_gene-microbe_mt_kegg_MICROBE.csv")
write.csv(cca_output$cors,'CSTIV_host_mt_gene-microbe_mt_kegg_CORR.csv')


