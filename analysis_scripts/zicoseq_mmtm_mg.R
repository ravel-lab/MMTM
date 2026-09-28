# Load GUniFrac library
library(GUniFrac)


#############################################################
# Script to compare the relative abundance of different taxa 
#  between the cervical and vaginal samples in the metagenomics
#############################################################

#load the metagenomics taxa abundance data determined by mapping reads to VIRGO2
mmtm_mg_comp <- read.csv("/Users/ichaudry/Documents/OneDrive - University of Maryland School of Medicine/research/ravel_lab/MMTM/analysis/00_mg_analysis/virgo2/taxa_filt_mmtm_mg_virgo2.csv")
rownames(mmtm_mg_comp) <- mmtm_mg_comp$Taxa
mmtm_mg_comp <- mmtm_mg_comp[,-1]
mmtm_mg_comp <- mmtm_mg_comp[rowSums(mmtm_mg_comp[]) > 0,]

#load the metadata for each sample
mmtm_metadata <- read.csv("/Users/ichaudry/Documents/OneDrive - University of Maryland School of Medicine/research/ravel_lab/MMTM/analysis/mmtm_metadata.csv")
rownames(mmtm_metadata) <- mmtm_metadata$IGS.LABEL
mmtm_metadata <- mmtm_metadata[,-1]


#load in the cst assingments per sample
sample_csts <- read.csv('/Users/ichaudry/Documents/OneDrive - University of Maryland School of Medicine/research/ravel_lab/MMTM/analysis/cst_analysis/mmtm_cst_assignments_taxa_corr.csv',sep=',',row.names = 1, header = TRUE)


#filter out the metadata dataframe to only include the samples that have composition data
samples = colnames(mmtm_mg_comp)
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
mmtm_mg_comp <- mmtm_mg_comp[, colnames(mmtm_mg_comp) %in% samples, drop=FALSE]



mmtm_mg_comp <- as.matrix(mmtm_mg_comp)
mmtm.ZicoSeq.obj <- ZicoSeq(meta.dat = mmtm_metadata, feature.dat = mmtm_mg_comp, 
                       grp.name = 'Cervix.Vagina', adj.name=c('CST','Ext.Participant.ID', 'Timepoint'),feature.dat.type = "proportion",
                       # Filter to remove rare taxa
                       prev.filter = 0.2, mean.abund.filter = 0.00008,  
                        min.prop = 0,
                       # Winsorization to replace outliers
                       is.winsor = TRUE, outlier.pct = 0.03, winsor.end = 'top',
                       # Posterior sampling 
                       is.post.sample = TRUE, post.sample.no = 25, 
                       # Use the square-root transformation
                       link.func = list(function (x) x^0.5), stats.combine.func = max,
                       # Permutation-based multiple testing correction
                       perm.no = 9999,  strata =  NULL, 
                       # Reference-based multiple stage normalization
                       ref.pct = 0.5, stage.no = 6, excl.pct = 0.2,
                       # Family-wise error rate control
                       is.fwer = TRUE, verbose = TRUE, return.feature.dat = TRUE)



ZicoSeq.plot(mmtm.ZicoSeq.obj, pvalue.type = 'p.adj.fdr', cutoff = .05, text.size = 10,
             out.dir = 'NULL', width = 10, height = 6)


ZicoSeq.plot(mmtm.ZicoSeq.obj, pvalue.type = 'p.adj.fdr', cutoff = .05, text.size = 10,
             out.dir = '/Users/ichaudry/Documents/OneDrive - University of Maryland School of Medicine/research/ravel_lab/MMTM/analysis/00_mg_analysis/zicoseq/mmtm_mg_taxa_cv_zicoseq_results', width = 10, height = 6)

coefs <- t(as.data.frame(mmtm.ZicoSeq.obj$coef.list[1]))
to_write <- as.data.frame(mmtm.ZicoSeq.obj[c('R2','p.adj.fdr')])
colnames(to_write)[1] <- 'R2'

to_write <- cbind(to_write, coefs)


write.csv(to_write, "/Users/ichaudry/Documents/OneDrive - University of Maryland School of Medicine/research/ravel_lab/MMTM/analysis/00_mg_analysis/zicoseq/mmtm_mg_taxa_cv_zicoseq_results.csv")


#############################################################
# Script to compare the relative abundance of different genes 
#  between the cervical and vaginal samples in the metagenomics
#############################################################

#load the metagenomics gene level abundance data determined by mapping reads to VIRGO2
mmtm_mg_comp <- read.csv("/Users/ichaudry/Documents/OneDrive - University of Maryland School of Medicine/research/ravel_lab/MMTM/analysis/00_mg_analysis/virgo2/taxa_filt_mmtm_mg_virgo2_NR_Taxa_103023.csv", header=TRUE, row.names=1)

#drop the category and length columns from the mapping table
mmtm_mg_comp <- subset(mmtm_mg_comp, select=-c(Taxa, Cat, Length))


#load in the cst assingments per sample
sample_csts <- read.csv('/Users/ichaudry/Documents/OneDrive - University of Maryland School of Medicine/research/ravel_lab/MMTM/analysis/cst_analysis/mmtm_cst_assignments_taxa_corr.csv',sep=',',row.names = 1, header = TRUE)


#load the metadata for each sample
mmtm_metadata <- read.csv("/Users/ichaudry/Documents/OneDrive - University of Maryland School of Medicine/research/ravel_lab/MMTM/analysis/mmtm_metadata.csv", sep=',', row.names = 1)

#filter out the metadata dataframe to only include the samples that have composition data
samples = colnames(mmtm_mg_comp)
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
mmtm_mg_comp <- mmtm_mg_comp[, colnames(mmtm_mg_comp) %in% samples, drop=FALSE]



#drop genes with all 0s
mmtm_mg_comp <- mmtm_mg_comp[rowSums(mmtm_mg_comp[]) > 0,]

mmtm_mg_comp <- as.matrix(mmtm_mg_comp)
mmtm.genes.ZicoSeq.obj <- ZicoSeq(meta.dat = mmtm_metadata, feature.dat = mmtm_mg_comp, 
                                  grp.name = 'Cervix.Vagina',adj.name = c('CST', 'Ext.Participant.ID'),feature.dat.type = "count",
                                  # Filter to remove rare taxa
                                  prev.filter = 0.25, mean.abund.filter = 0,  
                                  min.prop = 0,
                                  # Winsorization to replace outliers
                                  is.winsor = TRUE, outlier.pct = 0.03, winsor.end = 'top',
                                  # Posterior sampling 
                                  is.post.sample = TRUE, post.sample.no = 25, 
                                  # Use the square-root transformation
                                  link.func = list(function (x) x^0.5), stats.combine.func = max,
                                  # Permutation-based multiple testing correction
                                  perm.no = 99,  strata = NULL, 
                                  # Reference-based multiple stage normalization
                                  ref.pct = 0.5, stage.no = 6, excl.pct = 0.2,
                                  # Family-wise error rate control
                                  is.fwer = TRUE, verbose = TRUE, return.feature.dat = TRUE)



ZicoSeq.plot(mmtm.genes.ZicoSeq.obj, pvalue.type = 'p.adj.fdr', cutoff = .05, text.size = 10,
             out.dir = NULL, width = 10, height = 6)



ZicoSeq.plot(mmtm.genes.ZicoSeq.obj, pvalue.type = 'p.adj.fdr', cutoff = .05, text.size = 10,
             out.dir = '/Users/ichaudry/Documents/OneDrive - University of Maryland School of Medicine/research/ravel_lab/MMTM/analysis/00_mg_analysis/zicoseq/mmtm_mg_gene_cv_zicoseq_results', width = 10, height = 6)


coefs <- t(as.data.frame(mmtm.genes.ZicoSeq.obj$coef.list[1]))
to_write <- as.data.frame(mmtm.genes.ZicoSeq.obj[c('R2','p.adj.fdr')])
colnames(to_write)[1] <- 'R2'

to_write <- cbind(to_write, coefs)
write.csv(to_write, "/Users/ichaudry/Documents/OneDrive - University of Maryland School of Medicine/research/ravel_lab/MMTM/analysis/00_mg_analysis/zicoseq/mmtm_mg_gene_cv_zicoseq_results.csv")


#write the results to a csv
#w <- data.frame(mmtm.genes.ZicoSeq.obj$p.raw)
#x <- data.frame(mmtm.genes.ZicoSeq.obj$p.adj.fdr)
#y <- data.frame(mmtm.genes.ZicoSeq.obj$R2)
#z <- data.frame(mmtm.genes.ZicoSeq.obj$coef.list)
#z <- t(z)

#to_write <- cbind(w,x,y,z)
#write.csv(to_write,"/Users/ichaudry/Documents/OneDrive - University of Maryland School of Medicine/research/ravel_lab/MMTM/analysis/00_mg_analysis/zicoseq/mmtm_mg_gene_zicoseq_results.csv")




#############################################################
# Script to compare the relative abundance of different KEGG categories 
#  between the cervical and vaginal samples in the metagenomics
#############################################################

#load the metagenomics KEGG abundance data determined by mapping reads to VIRGO2
mmtm_mg_comp <- read.csv("/Users/ichaudry/Documents/OneDrive - University of Maryland School of Medicine/research/ravel_lab/MMTM/analysis/00_mg_analysis/virgo2/taxa_filt_mmtm_mg_virgo2_kegg.csv")
rownames(mmtm_mg_comp) <- mmtm_mg_comp$KEGG
mmtm_mg_comp <- mmtm_mg_comp[,-1]


#load the metadata for each sample
mmtm_metadata <- read.csv("/Users/ichaudry/Documents/OneDrive - University of Maryland School of Medicine/research/ravel_lab/MMTM/analysis/mmtm_metadata.csv")
rownames(mmtm_metadata) <- mmtm_metadata$IGS.LABEL
mmtm_metadata <- mmtm_metadata[,-1]


#load in the cst assingments per sample
sample_csts <- read.csv('/Users/ichaudry/Documents/OneDrive - University of Maryland School of Medicine/research/ravel_lab/MMTM/analysis/cst_analysis/mmtm_cst_assignments_taxa_corr.csv',sep=',',row.names = 1, header = TRUE)


#filter out the metadata dataframe to only include the samples that have composition data
samples = colnames(mmtm_mg_comp)
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
mmtm_mg_comp <- mmtm_mg_comp[, colnames(mmtm_mg_comp) %in% samples, drop=FALSE]


mmtm_mg_comp <- mmtm_mg_comp[rowSums(mmtm_mg_comp[]) > 0,]
mmtm_mg_comp <- as.matrix(mmtm_mg_comp)
mmtm.ZicoSeq.obj <- ZicoSeq(meta.dat = mmtm_metadata, feature.dat = mmtm_mg_comp, 
                            grp.name = 'Cervix.Vagina', adj.name=c('CST','Ext.Participant.ID'),feature.dat.type = "proportion",
                            # Filter to remove rare taxa
                            prev.filter = 0.2, mean.abund.filter = 0,  
                            min.prop = 0,
                            # Winsorization to replace outliers
                            is.winsor = TRUE, outlier.pct = 0.03, winsor.end = 'top',
                            # Posterior sampling 
                            is.post.sample = TRUE, post.sample.no = 25, 
                            # Use the square-root transformation
                            link.func = list(function (x) x^0.5), stats.combine.func = max,
                            # Permutation-based multiple testing correction
                            perm.no = 999,  strata =  NULL, 
                            # Reference-based multiple stage normalization
                            ref.pct = 0.5, stage.no = 6, excl.pct = 0.2,
                            # Family-wise error rate control
                            is.fwer = TRUE, verbose = TRUE, return.feature.dat = TRUE)



ZicoSeq.plot(mmtm.ZicoSeq.obj, pvalue.type = 'p.adj.fdr', cutoff = .05, text.size = 10,
             out.dir = 'NULL', width = 10, height = 6)


ZicoSeq.plot(mmtm.ZicoSeq.obj, pvalue.type = 'p.adj.fdr', cutoff = .05, text.size = 10,
             out.dir = '/Users/ichaudry/Documents/OneDrive - University of Maryland School of Medicine/research/ravel_lab/MMTM/analysis/00_mg_analysis/zicoseq/mmtm_mg_kegg_cv_zicoseq_results', width = 10, height = 6)

coefs <- t(as.data.frame(mmtm.ZicoSeq.obj$coef.list[1]))
to_write <- as.data.frame(mmtm.ZicoSeq.obj[c('R2','p.adj.fdr')])
colnames(to_write)[1] <- 'R2'

to_write <- cbind(to_write, coefs)

write.csv(to_write, "/Users/ichaudry/Documents/OneDrive - University of Maryland School of Medicine/research/ravel_lab/MMTM/analysis/00_mg_analysis/zicoseq/mmtm_mg_kegg_cv_zicoseq_results.csv")



#################################################################
#################################################################
######### BELOW THIS IS ALL OVER TIME COMPARISONS ###############
#################################################################
#################################################################

#############################################################
# Script to compare the relative abundance of different taxa 
#  between TIMEPOINTS in the cervix.
#############################################################

#load the metagenomics taxa abundance data determined by mapping reads to VIRGO2
mmtm_mg_comp <- read.csv("/Users/ichaudry/Documents/OneDrive - University of Maryland School of Medicine/research/ravel_lab/MMTM/analysis/00_mg_analysis/virgo2/taxa_filt_mmtm_mg_virgo2.csv")
rownames(mmtm_mg_comp) <- mmtm_mg_comp$Taxa
mmtm_mg_comp <- mmtm_mg_comp[,-1]


#load the metadata for each sample
mmtm_metadata <- read.csv("/Users/ichaudry/Documents/OneDrive - University of Maryland School of Medicine/research/ravel_lab/MMTM/analysis/mmtm_metadata.csv")
rownames(mmtm_metadata) <- mmtm_metadata$IGS.LABEL
mmtm_metadata <- mmtm_metadata[,-1]


#load in the cst assingments per sample
sample_csts <- read.csv('/Users/ichaudry/Documents/OneDrive - University of Maryland School of Medicine/research/ravel_lab/MMTM/analysis/cst_analysis/mmtm_cst_assignments_taxa_corr.csv',sep=',',row.names = 1, header = TRUE)


#filter out the metadata dataframe to only include the samples that have composition data
samples = colnames(mmtm_mg_comp)
mmtm_metadata <- mmtm_metadata[row.names(mmtm_metadata) %in% samples, , drop=FALSE]
mmtm_metadata$Cervix.Vagina <- as.factor(mmtm_metadata$Cervix.Vagina)
mmtm_metadata$Timepoint <- as.factor(mmtm_metadata$Timepoint)
mmtm_metadata$Ext.Participant.ID <- as.factor(mmtm_metadata$Ext.Participant.ID)

#merge the CST information with the metadata
mmtm_metadata <- merge(sample_csts, mmtm_metadata, by='row.names' )
rownames(mmtm_metadata) <- mmtm_metadata$Row.names
mmtm_metadata$CST <- as.factor(mmtm_metadata$CST)

#filter for only 1 location samples
mmtm_metadata <- mmtm_metadata[mmtm_metadata$Cervix.Vagina == 'V',]

#filter out the control samples from the composition data by only keeping the samples that also have metadata associated with them
samples = rownames(mmtm_metadata)
mmtm_mg_comp <- mmtm_mg_comp[, colnames(mmtm_mg_comp) %in% samples, drop=FALSE]


mmtm_mg_comp <- mmtm_mg_comp[rowSums(mmtm_mg_comp[]) > 0,]
mmtm_mg_comp <- as.matrix(mmtm_mg_comp)
mmtm.ZicoSeq.obj <- ZicoSeq(meta.dat = mmtm_metadata, feature.dat = mmtm_mg_comp, 
                            grp.name = 'Timepoint', adj.name=c('CST'),feature.dat.type = "proportion",
                            # Filter to remove rare taxa
                            prev.filter = 0, mean.abund.filter = 0.00001,  
                            min.prop = 0,
                            # Winsorization to replace outliers
                            is.winsor = TRUE, outlier.pct = 0.03, winsor.end = 'top',
                            # Posterior sampling 
                            is.post.sample = TRUE, post.sample.no = 25, 
                            # Use the square-root transformation
                            link.func = list(function (x) x^0.5), stats.combine.func = max,
                            # Permutation-based multiple testing correction
                            perm.no = 20000,  strata =  NULL, 
                            # Reference-based multiple stage normalization
                            ref.pct = 0.5, stage.no = 6, excl.pct = 0.2,
                            # Family-wise error rate control
                            is.fwer = TRUE, verbose = TRUE, return.feature.dat = TRUE)



ZicoSeq.plot(mmtm.ZicoSeq.obj, pvalue.type = 'p.adj.fdr', cutoff = .05, text.size = 10,
             out.dir = 'NULL', width = 10, height = 6)


ZicoSeq.plot(mmtm.ZicoSeq.obj, pvalue.type = 'p.adj.fdr', cutoff = .05, text.size = 10,
             out.dir = '/Users/ichaudry/Documents/OneDrive - University of Maryland School of Medicine/research/ravel_lab/MMTM/analysis/00_mg_analysis/zicoseq/mmtm_mg_taxa_tp_zicoseq_results', width = 10, height = 6)

coefs <- t(as.data.frame(mmtm.ZicoSeq.obj$coef.list[1]))
to_write <- as.data.frame(mmtm.ZicoSeq.obj[c('R2','p.adj.fdr')])
colnames(to_write)[1] <- 'R2'

to_write <- cbind(to_write, coefs)


write.csv(to_write, "/Users/ichaudry/Documents/OneDrive - University of Maryland School of Medicine/research/ravel_lab/MMTM/analysis/00_mg_analysis/zicoseq/mmtm_mg_taxa_tp_zicoseq_results.csv")




#############################################################
# Script to compare the relative abundance of different genes 
#  between TIMEPOINTS
#############################################################

#load the metagenomics gene level abundance data determined by mapping reads to VIRGO2
mmtm_mg_comp <- read.csv("/Users/ichaudry/Documents/OneDrive - University of Maryland School of Medicine/research/ravel_lab/MMTM/analysis/00_mg_analysis/virgo2/taxa_filt_mmtm_mg_virgo2_NR_Taxa_103023.csv", header=TRUE, row.names=1)

#drop the category and length columns from the mapping table
mmtm_mg_comp <- subset(mmtm_mg_comp, select=-c(Taxa, Cat, Length))


#load in the cst assingments per sample
sample_csts <- read.csv('/Users/ichaudry/Documents/OneDrive - University of Maryland School of Medicine/research/ravel_lab/MMTM/analysis/cst_analysis/mmtm_cst_assignments_taxa_corr.csv',sep=',',row.names = 1, header = TRUE)


#load the metadata for each sample
mmtm_metadata <- read.csv("/Users/ichaudry/Documents/OneDrive - University of Maryland School of Medicine/research/ravel_lab/MMTM/analysis/mmtm_metadata.csv", sep=',', row.names = 1)

#filter out the metadata dataframe to only include the samples that have composition data
samples = colnames(mmtm_mg_comp)
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
mmtm_mg_comp <- mmtm_mg_comp[, colnames(mmtm_mg_comp) %in% samples, drop=FALSE]



#drop genes with all 0s
mmtm_mg_comp <- mmtm_mg_comp[rowSums(mmtm_mg_comp[]) > 0,]

mmtm_mg_comp <- as.matrix(mmtm_mg_comp)
mmtm.genes.ZicoSeq.obj <- ZicoSeq(meta.dat = mmtm_metadata, feature.dat = mmtm_mg_comp, 
                                  grp.name = 'Timepoint',adj.name = c('CST', 'Ext.Participant.ID'),feature.dat.type = "count",
                                  # Filter to remove rare taxa
                                  prev.filter = 0.25, mean.abund.filter = 0.00001,  
                                  min.prop = 0,
                                  # Winsorization to replace outliers
                                  is.winsor = TRUE, outlier.pct = 0.03, winsor.end = 'top',
                                  # Posterior sampling 
                                  is.post.sample = TRUE, post.sample.no = 25, 
                                  # Use the square-root transformation
                                  link.func = list(function (x) x^0.5), stats.combine.func = max,
                                  # Permutation-based multiple testing correction
                                  perm.no = 99,  strata = NULL, 
                                  # Reference-based multiple stage normalization
                                  ref.pct = 0.5, stage.no = 6, excl.pct = 0.2,
                                  # Family-wise error rate control
                                  is.fwer = TRUE, verbose = TRUE, return.feature.dat = TRUE)



ZicoSeq.plot(mmtm.genes.ZicoSeq.obj, pvalue.type = 'p.adj.fdr', cutoff = .05, text.size = 10,
             out.dir = NULL, width = 10, height = 6)



ZicoSeq.plot(mmtm.genes.ZicoSeq.obj, pvalue.type = 'p.adj.fdr', cutoff = .05, text.size = 10,
             out.dir = '/Users/ichaudry/Documents/OneDrive - University of Maryland School of Medicine/research/ravel_lab/MMTM/analysis/00_mg_analysis/zicoseq/mmtm_mg_gene_tp_zicoseq_results', width = 10, height = 6)




coefs <- t(as.data.frame(mmtm.genes.ZicoSeq.obj$coef.list[1]))
to_write <- as.data.frame(mmtm.genes.ZicoSeq.obj[c('R2','p.adj.fdr')])
colnames(to_write)[1] <- 'R2'

to_write <- cbind(to_write, coefs)
write.csv(to_write, "/Users/ichaudry/Documents/OneDrive - University of Maryland School of Medicine/research/ravel_lab/MMTM/analysis/00_mg_analysis/zicoseq/mmtm_mg_genes_tp_zicoseq_results.csv")



#write the results to a csv
#w <- data.frame(mmtm.genes.ZicoSeq.obj$p.raw)
#x <- data.frame(mmtm.genes.ZicoSeq.obj$p.adj.fdr)
#y <- data.frame(mmtm.genes.ZicoSeq.obj$R2)
#z <- data.frame(mmtm.genes.ZicoSeq.obj$coef.list)
#z <- t(z)

#to_write <- cbind(w,x,y,z)
#write.csv(to_write,"/Users/ichaudry/Documents/OneDrive - University of Maryland School of Medicine/research/ravel_lab/MMTM/analysis/00_mg_analysis/zicoseq/mmtm_mg_gene_zicoseq_results.csv")


#############################################################
# Script to compare the relative abundance of different KEGG categories 
#  between TIMEPOINTS
#############################################################

#load the metagenomics KEGG abundance data determined by mapping reads to VIRGO2
mmtm_mg_comp <- read.csv("/Users/ichaudry/Documents/OneDrive - University of Maryland School of Medicine/research/ravel_lab/MMTM/analysis/00_mg_analysis/virgo2/taxa_filt_mmtm_mg_virgo2_kegg.csv")
rownames(mmtm_mg_comp) <- mmtm_mg_comp$KEGG
mmtm_mg_comp <- mmtm_mg_comp[,-1]


#load the metadata for each sample
mmtm_metadata <- read.csv("/Users/ichaudry/Documents/OneDrive - University of Maryland School of Medicine/research/ravel_lab/MMTM/analysis/mmtm_metadata.csv")
rownames(mmtm_metadata) <- mmtm_metadata$IGS.LABEL
mmtm_metadata <- mmtm_metadata[,-1]


#load in the cst assingments per sample
sample_csts <- read.csv('/Users/ichaudry/Documents/OneDrive - University of Maryland School of Medicine/research/ravel_lab/MMTM/analysis/cst_analysis/mmtm_cst_assignments_taxa_corr.csv',sep=',',row.names = 1, header = TRUE)


#filter out the metadata dataframe to only include the samples that have composition data
samples = colnames(mmtm_mg_comp)
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
mmtm_mg_comp <- mmtm_mg_comp[, colnames(mmtm_mg_comp) %in% samples, drop=FALSE]


mmtm_mg_comp <- mmtm_mg_comp[rowSums(mmtm_mg_comp[]) > 0,]
mmtm_mg_comp <- as.matrix(mmtm_mg_comp)
mmtm.ZicoSeq.obj <- ZicoSeq(meta.dat = mmtm_metadata, feature.dat = mmtm_mg_comp, 
                            grp.name = 'Timepoint', adj.name=c('CST','Ext.Participant.ID'),feature.dat.type = "proportion",
                            # Filter to remove rare taxa
                            prev.filter = 0.2, mean.abund.filter = 0,  
                            min.prop = 0,
                            # Winsorization to replace outliers
                            is.winsor = TRUE, outlier.pct = 0.03, winsor.end = 'top',
                            # Posterior sampling 
                            is.post.sample = TRUE, post.sample.no = 25, 
                            # Use the square-root transformation
                            link.func = list(function (x) x^0.5), stats.combine.func = max,
                            # Permutation-based multiple testing correction
                            perm.no = 999,  strata =  NULL, 
                            # Reference-based multiple stage normalization
                            ref.pct = 0.5, stage.no = 6, excl.pct = 0.2,
                            # Family-wise error rate control
                            is.fwer = TRUE, verbose = TRUE, return.feature.dat = TRUE)



ZicoSeq.plot(mmtm.ZicoSeq.obj, pvalue.type = 'p.adj.fdr', cutoff = .05, text.size = 10,
             out.dir = 'NULL', width = 10, height = 6)


ZicoSeq.plot(mmtm.ZicoSeq.obj, pvalue.type = 'p.adj.fdr', cutoff = .05, text.size = 10,
             out.dir = '/Users/ichaudry/Documents/OneDrive - University of Maryland School of Medicine/research/ravel_lab/MMTM/analysis/00_mg_analysis/zicoseq/mmtm_mg_kegg_tp_zicoseq_results', width = 10, height = 6)

coefs <- t(as.data.frame(mmtm.ZicoSeq.obj$coef.list[1]))
to_write <- as.data.frame(mmtm.ZicoSeq.obj[c('R2','p.adj.fdr')])
colnames(to_write)[1] <- 'R2'

to_write <- cbind(to_write, coefs)


write.csv(to_write, "/Users/ichaudry/Documents/OneDrive - University of Maryland School of Medicine/research/ravel_lab/MMTM/analysis/00_mg_analysis/zicoseq/mmtm_mg_kegg_tp_zicoseq_results.csv")















