rm(list=ls())
library(data.table)
library(plyr)
library(ggplot2)
library(viridis)
library(cowplot)
library(stringr)
library(splines)
setwd('set-path-appropriately')
source("./registration_functions.R")

# R version 4.0.0 (2020-04-24)
# Platform: x86_64-apple-darwin17.0 (64-bit)
# Running under: macOS Catalina 10.15.7

# data.table_1.13.2
# plyr_1.8.6
# ggplot2_3.3.2 
# viridis_0.5.1
# cowplot_1.0.0
# stringr_1.4.0
# splines_4.0.0


##### SET PARAMETERS:

# parameters for performing registration:
do.initial.rescale <- 'false'  # 'rescale' to scale gene expression prior to registration.
do.register.rescale <- 'rescale'  # 'rescale' to scale gene expression using only overlapping timepoints
                                  # points during registration.
stretch = c(2, 1.5, 1) # candidate registration stretch factors to apply to Col-0
min.num.overlapping.points = 4 # will only consider shifts which leave at least this many overlapping timepoints 
                              # after applying the registration function.
shift.extreme  = 4 # the absolute maximum value which can be applied as a shift to gene expression timecourse (days).
transformed.timecourse = 'Col0' # the name of the timecourse to apply the registration function to.
num.shuffled = 1000  # the number of random permutations to use to test real results for significance.
shuffle.type <- 'shuffle.genes' # method for random permutation of the data to assess significance of real results. 
                                # 'shuffle.genes' : randomly reallocate orthologues in Arabdidopsis compared to 
                                # each R-o-18 gene.
in.data.path <- './data_for_registration.rds'  # path to gene expression data
out.dir <- './registration_out/'  # path to directory to save results

##### END PARAMETERS


# create out.dir if doesn't exist
dir.create(out.dir, showWarnings = F)

# convert flags for rescaling to bools.
if (do.initial.rescale=='rescale') {
  initial.rescale <- TRUE
} else {
  initial.rescale <- FALSE
}
if (do.register.rescale=='rescale') {
  should.rescale <- TRUE
} else {
  should.rescale <- FALSE
}


##### LOAD GENE EXPRESSION DATA:
# load expression data for compared genes.

# mean.df : the mean gene expression of each gene in each genotype at each timepoint. 
# Columns:
# "locus_name" : brassica gene identifier. (identify of brassica genes. identity of brassica orthologue of Arabidopsis genes)
# "accession" : genotype - "Col0" or "Ro18"
# "timepoint" : the days post germination the sample taken
# "mean.cpm" : mean of TMM normalised counts per million amoung biological replicates.

# all.data.df : all replicates of gene expression in each genotype at each timepoint.
# Columns: 
# "locus_name" : brassica gene identifier. (identify of brassica genes. identity of brassica orthologue of Arabidopsis genes)
# "accession" : genotype - "Col0" or "Ro18"
# "timepoint" : the days post germination the sample taken
# "group" : library identifier, format "genotype-timepoint-replicate"
# "mean.cpm" : TMM normalised counts per million in each biological replicate (group).

L <- readRDS(in.data.path) # sumBrassica copy flag is within load_mean.df()
mean.df <- L[[1]] # will compare Col0 vs Ro18 based on accession column
all.data.df <- L[[2]]


################ cut down to few genes for testing #####################
# test.genes <- unique(mean.df$locus_name)[1:50]
# mean.df <- mean.df[mean.df$locus_name %in% test.genes,]
# all.data.df <- all.data.df[all.data.df$locus_name %in% test.genes, ]
########################################################################





# record original, unpermuted data
unshuffled.data <- copy(mean.df)
unshuffled.all.data <- copy(all.data.df)

# iterate over real, and randomly permuted data for registration and distance calculations.
real.and.shuffled <- c('real', 1:num.shuffled)
for (i in real.and.shuffled) {

  ##### RANDOMLY PERMUTE THE DATA IF DOING SHUFFLING FOR STATISTICAL SIGNIFICANCE TESTING OF REAL REGISTRATION RESULT
  if (i == 'real') {
    mean.df <- copy(unshuffled.data)
    all.data.df <- copy(unshuffled.all.data)
  } else {
    L <- shuffle_ro18_gene_names(mean.df, all.data.df)
    mean.df <- L[[1]]
    all.data.df <- L[[2]]
  }
   
  
  ##### REGISTER THE GENE EXPRESSION DATA
  O <- prepare_scaled_and_registered_data(mean.df, all.data.df, stretch=stretch, initial.rescale, should.rescale, min.num.overlapping.points, 
                                          shift.extreme, transformed.timecourse)

  mean.df <- O[['mean.df']] # mean.df is unchanged by "prepare_scaled_and_registered_data()
  mean.df.sc <- O[['mean.df.sc']] # mean.df.sc : Identical to mean.df, with additional column "sc.mean.cpm" : gene 
                                  # expression data scaled using "scale(center=T, scale=T)"
  imputed.mean.df <- O[['imputed.mean.df']] # imputed.mean.df is registered expression data. After registration function was applied to Col0 
                                            # timepoints, piecewise linear regression carried out to extimate expression at same timepoints 
                                            # have Ro18 observations. 
                                            # "shifted time" column is timepoint after registration applied.
                                            # "is.registered" boolean column of whether model selection found 
                                            # registration model (TRUE) or no-registration model (FALSE) is more supported.
  all.shifts <- O[['all.shifts']] # table of candidate registraions applied, and score for each.
                                  # "gene": brassica gene id
                                  # "stretch : stretch factor applied to Col0 timecourse in registration function
                                  # "shift" : shift applied to Col0 timecourse (days) in registration function
                                  # "score" : mean squared difference between Col0 and Ro18 after registration function applied
  model.comparison <- O[['model.comparison']] # table comparing the optimal registration function for each gene (based on
                                              # minimizing all.shifts$score) to model with no registration applied.
                                              # "gene" : brassica gene id.
                                              # "seperate.AIC" : Akaike Information Score for spline modelling 
                                              # of gene expression in Ro18 and Sarisha14 seperately.
                                              # "registered.AIC" : Akaike Information Score for single model of gene expression 
                                              # in Ro18 and registered Col0 expression.
                                              # "seperate.BIC" : Bayesian Information Criterion score for spline modelling 
                                              # of gene expression in Ro18 and Sarisha14 seperately.
                                              # "registered.BIC" : Bayesian Information Criterion score for ingle model of gene expression 
                                              # in Ro18 and registered Col0 expression.
                                              # "stretch" : stretch factor applied to Col0 timecourse in optimal registration function
                                              # "shift" : shift applied to Col0 timecourse (days) in optimal registration function
                                              # "BIC.registered.is.better" : boolean for whether regitration model is more supported
                                              # than seperate model by BIC criterion.
                                              # "AIC.registered.is.better" : boolean for whether regitration model is more supported
                                              # than seperate model by AIC criterion.
                                              # "ABIC.registered.is.better" : boolean for whether regitration model is more supported
                                              # than seperate model by both BIC & AIC criterion.
                                            
  saveRDS(O, paste0(out.dir, 'registration_result_', i, '.rds'))
  
  
  ##### CALCULATE DISTANCES BETWEEN THE SAMPLES ####
  # calculate distance between samples based on gene expression. 
  # Mean expression amoung biological replicates is used to calculate distance. 
  # Distance metric is mean euclidean distance between expression of orthologous gene pairs. 
  O <- calculate_between_sample_distance(mean.df, mean.df.sc, imputed.mean.df)
  
  D.mean <- O[['D.mean']]  # distance using gene expression averaged between replicates.
  D.scaled <- O[['D.scaled']]  # distance using scaled gene expression
  D.registered <- O[['D.registered']]  # distance after registration
  
  saveRDS(O, paste0(out.dir, 'distance_result_', i, '.rds'))
}

