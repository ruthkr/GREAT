rm(list=ls())
library(data.table)
library(ggplot2)
library(stringr)
library(edgeR)
library(GO.db) # GO annotation package
library(org.At.tair.db) # arabidopsis annotation package
library(clusterProfiler)

sessionInfo()
# R version 4.0.0 (2020-04-24)
# Platform: x86_64-apple-darwin17.0 (64-bit)
# Running under: macOS Catalina 10.15.

# ggplot2_3.3.2
# data.table_1.13.2
# stringr_1.4.0
# edgeR_3.30.3
# clusterProfiler_3.16.1
# org.At.tair.db_3.11.4
# GO.db_3.11.4
# AnnotationDbi_1.50.3


##### READ AND PREPARE DATA

# df is data.table object with columns:
# CDS.model : gene identifiers
# sample_id : sample / library unique identifiers
# accession : "sarisha14" or "Ro18" identifier
# timepoint : days post germination sample taken 
# est_counts : Stringtie estimated gene expression counts.
df <- load_data()

# cut down to timepoints for pairwise comparison
early.df <- df[(df$accession=='sarisha14' & df$timepoint==9) |
                 df$accession=='Ro18' & df$timepoint==11,]
early.df <- unique(early.df[, c('CDS.model', 'sample_id', 'est_counts')])
early.df.c <- dcast(early.df, CDS.model~sample_id, value.var='est_counts')

# if genes detected in one genotype, but not the other, count is NA. Replace with 0
early.df.c[is.na(early.df.c)] <- 0



##### DIFF. EXPRESSION ANALYSIS, 
# following https://bioinformatics-core-shared-training.github.io/cruk-bioinf-sschool/Day3/rnaSeq_DE.pdf

# make into EdgeR object
early.df.m <- as.matrix(early.df.c[,2:ncol(early.df.c)])
rownames(early.df.m) <- early.df.c$CDS.model
dgList <- DGEList(counts= early.df.m, genes=rownames(early.df.m))

# filtering - retain genes which at least 1cpm in at least 2 samples
countsPerMillion <- cpm(dgList)
countCheck <- countsPerMillion > 1
keep <- which(rowSums(countCheck) >= 2)
dgList <- dgList[keep,]
    
# between library normalisation        
dgList <- calcNormFactors(dgList, method='TMM')

# setup experimental design to use for model fitting
sampleType <- rep('ro18', ncol(dgList))
sampleType[grep('sari', colnames(dgList))] <- 'sari14'
sampleType <- as.factor(sampleType)
designMat <- model.matrix(~sampleType)

# Estimating dispersions
dgList <- estimateGLMCommonDisp(dgList, design=designMat)
dgList <- estimateGLMTrendedDisp(dgList, design=designMat)
dgList <- estimateGLMTagwiseDisp(dgList, design=designMat)

# Differential expression
fit <- glmFit(dgList, designMat)
lrt <- glmLRT(fit, coef='sampleTypesari14')




##### GO ENRICHMENT ANALYSIS
# reformat glmLRT output
DE.results <- lrt$table
DE.results$gene <- rownames(DE.results)
deGenes <- decideTestsDGE(lrt, adjust.method='BH', p=0.05) # significantly DE genes. 
                                                           # Benjamini-Hochberg correction for multiple testing.
                                                           # significance threshold 0.05.
deGenes.df <- data.frame('gene'=rownames(deGenes@.Data), 
                         'is.de'=deGenes@.Data)
DE.results <- merge(DE.results, deGenes.df, by='gene')
names(DE.results)[names(DE.results)=='sampleTypesari14'] <- 'diff.expressed.in.sari'

# map the diff. expressed Brassica genes to the Arabidopsis orthologoues TAIR ids
# (will use Arabidopsis GO info for enrichment test).
id.table <- load_orthologue_table() # table of brassica gene id : orthologous arabidopsis TAIR id
DE.results <- merge(DE.results, id.merged, by.x='gene', by.y='merged.id',
                    all.x=TRUE)

# do the GO enrichment test
background.genes <- DE.results$locus_name # vector of TAIR ids of arabidopsis orthologues of 
                                          # genes which could have been detected as differently expressed
                                          # because id detected in brassica, and because
                                          # orthologue in Arabidopsis is identified.
test.genes <- DE.results$locus_name[DE.results$diff.expressed.in.sari!=0]  # vector of TAIR ids of genes 
                                          # which ARE differently expressed.
diff_ego<- enrichGO(test.genes, 
                       org.At.tair.db, 
                       keyType='TAIR', 
                       ont='BP', 
                       universe=background.genes, 
                       minGSSize=3, # the minimum number of genes with the associated ontology in the universe list, for GO to be considered for possible enrichment
                       maxGSSize=500) # ditto for maximum
diff_GO <- diff_ego@result[,]


