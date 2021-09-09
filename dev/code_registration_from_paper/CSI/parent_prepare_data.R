rm(list=ls())

library(methods)
library(data.table)
library(kernlab)
library(ggplot2)
library(plyr)
library(cowplot)
library(viridis)

equality_test <- function(x1, x2, y1, y2, nStdDev=2) {
  # uses "Benavoli & Manglili" Gauss. Process method for pairwise equality tests.
  # is are the same, returns 1, if are not the same, returns 0.
  # nStdDev specifies how many standard deviations away from mean extimate of difference can include the 0 vector.
  
  # normalise observed data - nb that kernlab seems to really struggle with sd estimates if data isn't normalised first...
  y1.n <- scale(y1)
  y2.n <- scale(y2)
  #y1.n <- y1
  #y2.n <- y2
  
  # train gp1 and gp2
  GP1 <- gausspr(x1, y1.n, kernel='rbfdot', variance.model=T)
  GP2 <- gausspr(x2, y2.n, kernel='rbfdot', variance.model=T)
  
  # predict values for x_star mean, and standard deviation under both models
  x_star <- c(x1,x2)
  
  u_star1 <- predict(GP1, x_star)
  K_star1 <- predict(GP1, x_star, type="variance")
  u_star2 <- predict(GP2, x_star)
  K_star2 <- predict(GP2, x_star, type="variance")
  
  # sanity plots
  ggplot(data.frame('x'=x_star, 'mean'=u_star1, 'UL'=u_star1+nStdDev*sqrt(K_star1), 'LL'=u_star1-nStdDev*sqrt(K_star1)))+
    geom_point(data=data.frame('xt'=x1, 'yt'=y1.n), aes(x=xt, y=yt))+
    geom_line(aes(x=x, y=mean))+
    geom_ribbon(aes(x=x_star, ymin=LL, ymax=UL), alpha=0.2)
  
  ggplot(data.frame('x'=x_star, 'mean'=u_star2, 'UL'=u_star2+nStdDev*sqrt(K_star2), 'LL'=u_star2-nStdDev*sqrt(K_star2)))+
    geom_point(data=data.frame('xt'=x2, 'yt'=y2.n), aes(x=xt, y=yt))+
    geom_line(aes(x=x, y=mean))+
    geom_ribbon(aes(x=x_star, ymin=LL, ymax=UL), alpha=0.2)
  
  # get differences
  delta_u <- u_star1 - u_star2
  delta_K <- sqrt(K_star1 + K_star2) # takes variances and gives standard deviation
  
  # evaluate whether [0] vector in the 95% CI. -> nb NOT just that one of them is 0, as e.g. wavy function 
  # will meander through 0 at some points...
  UL <- delta_u + nStdDev*(delta_K)
  LL <- delta_u - nStdDev*(delta_K)
  
  # sanity plot of difference
  # library(ggplot2)
  # D <- data.frame('x'=x_star, 'delta_u'=delta_u, 'UL'=UL, 'LL'=LL)
  # p <- ggplot(D, aes(x=x, y=delta_u))+
  #   geom_line()+
  #   geom_ribbon(aes(ymin=LL, ymax=UL), alpha=0.2)
  # p
  
  
  # sanity plots of how G1, G2, difference look, and whether there's a difference
  plot.d <- data.frame
  
  null_vec <- rep(0,length(UL))
  
  # "outliers" is if at any timepoint, there is evidence that there is a difference between the 2 functions...
  outliers <- any(LL > 0 | UL < 0)
  #!outliers
  return( !outliers )
}

get_gene_expression = function(gene_list, expression_table) {
  # nb creates duplicate rows, when more than one symbol or description for each CDS.model
  #create empty df with the correct cols
  
  # if have some genes of interest
  if (!is.null(gene_list)) {
    list_of_dfs <- vector(mode='list', length=length(gene_list))
    for (i in 1:length(gene_list)) {
      id <- gene_list[i]
      lookup_column <- get_lookup_column(id)
      
      RESULTS = expression_table[ expression_table[[lookup_column]] == id, ]
      list_of_dfs[[i]] <- RESULTS
      
      if(nrow(RESULTS)==0) {
        print(paste("Sorry, I couldn't find any results for :", id))
      }
      
    }
    gene_expression <- rbind.fill(list_of_dfs)
    
  } else {
    # if no genes of interest, return whole table
    gene_expression <- expression_table
  }
  
  return(gene_expression)
}

get_lookup_column=function(id) {
  # take identifier, and return what column should try and look it up in
  
  id <- toupper(id)
  lookup_column <- NULL
  #work out what sort of id provided, and use to change the column searched for
  if(is_TAIR_gm(id)) {
    lookup_column <- "gene_id"
  } else if (is_TAIR_loc(id)) {
    lookup_column <- "locus_name"
  } else if (is_bancroft_id(id)) {
    lookup_column <- "CDS.model"
  } else {
    lookup_column <- "symbol"
  }
  print(paste0(id, " : lookup_cds_model... ...thinks id is: ", lookup_column))
  return(lookup_column)
}

is_TAIR_gm=function(id) {
  #test if string corresponds to TAIR AGI format for a gene model (i.e. .x), nb returns false for LOCI (i.e. without .X in it), 
  #instead requires a GENE model.
  return(all(grepl("AT[0-5M]G[0-9]+\\.[0-9]$", id)))
}

is_TAIR_loc=function(id) {
  #test if string corresponds to TAIR AGI format for a gene locus (i.e. NO .x), 
  return(all(grepl("AT[0-5M]G[0-9]+$", id)))
}

is_bancroft_id=function(id) {
  #test if string corresponds to Bancroft group assigned ID
  # e.g. Cab020598.1, BnaA01g05200D, Bo1g152060.1
  if (grepl("^CAB[0-9]+\\.[0-9]",id)) {
    return(TRUE)
  } else if (grepl("^BNA[AC][0-9]+G[0-9]",id)) {
    return(TRUE)
  } else if (grepl("^BO[0-9][r]?G[0-9]+",id)) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}

get_best_shift <- function(curr_sym, test) {
  min_shift <- 5
  
  # cut to get a single symbol
  test <- test[test$label==curr_sym, ]
  
  # get the mean expression vectors for the current gene
  AtVec <- test$mean.cpm[test$accession=='Col0']
  BrVec <- test$mean.cpm[test$accession=='Ro18']
  
  # get the arabidopsis, brassica timepoints will to compare to Col0 to
  # braMin = 1
  # atMin = 1
  # length = 6
  # 
  # cur_BrVec <- BrVec[seq(braMin,braMin+length)]
  # cur_AtVec <- AtVec[seq(atMin,atMin+length)]
  
  #data <- data.frame(at=cur_AtVec, br=cur_BrVec)
  
  # try sliding arabidopsis from left to right across brassica. Normalise expression by the considered timepoints (in each) mean
  # score inverse of mean distance between points - nb THIS STRETCHES, S.T. POINTS ARE LINED UP, EVEN IF WERE SAMPLED AT DIFFERENT TIMES!!!
  scores <- c()
  num_points <- 10
  for (num_points in seq(min_shift, length(BrVec))) {
    if (num_points <= length(AtVec)) {
      cur_at <- AtVec[seq(length(AtVec) - (num_points-1), length(AtVec))]
      cur_bra <- BrVec[seq(1, num_points)]
    } else {
      cur_at <- AtVec
      cur_bra <- BrVec[seq(1+num_points-length(AtVec), num_points)]
    }
    
    # rescale both
    cur_at <- cur_at / mean(cur_at)
    cur_bra <- cur_bra / mean(cur_bra)
    
    score <- 1 / mean((cur_at - cur_bra)^2)
    scores <- c(scores, score)
    
    adf <- data.frame(expression=cur_at,accession='at',time=seq(1,min(num_points, length(AtVec))))
    bdf <- data.frame(expression=cur_bra,accession='br',time=seq(1,min(num_points, length(AtVec))))
    df <- rbind(adf, bdf)
    ggplot(df, aes(x=time, y=expression, color=accession))+
      geom_point()+
      ggtitle(paste0(num_points, '->', score))
  }
  scores[is.na(scores)] <- 0 # if all 0, gives Na
  #scores
  best_ind <- which(scores==max(scores))
  best_num_points <- seq(min_shift, length(BrVec))[best_ind]
  
  return(best_num_points)
}

get_shifted_expression <- function(shift_results, exp) {
  group.bits <- tstrsplit(exp$group, '-')
  exp$accession <- group.bits[[1]]
  exp$tissue <- group.bits[[2]]
  exp$timepoint <- as.numeric(group.bits[[3]])
  
  cur_gene <- unique(shift_results$symbol)[2]
  shifted_exp <- list()
  for (cur_gene in unique(shift_results$symbol)) {
    
    # cut to get a single symbol
    test <- exp[exp$label==cur_gene, ]
    
    # get the expression vectors for the current gene
    atdf <- test[test$accession=='Col0']
    brdf <- test[test$accession=='Ro18']
    
    #AtVec <- atdf$norm.cpm #test$norm.cpm[test$accession=='Col0']
    #BrVec <- brdf$norm.cpm #test$norm.cpm[test$accession=='Ro18']
    
    num_points <- shift_results$num.points[shift_results$symbol==cur_gene]
    
    # get the timepoints which are considered to overlap. Will be the last "num_points" arabidopsis, and the first "num_points" brassica. 
    # if > 10, will be ALL ara timepoints!
    ara.timepoints <- sort(unique(test$timepoint[test$accession=='Col0']), decreasing=TRUE)
    if (num_points <= length(ara.timepoints)) {
      ara.time.overlap <- ara.timepoints[1:num_points]
    } else {
      ara.time.overlap <- ara.timepoints
    }
    bra.timepoints <- sort(unique(test$timepoint[test$accession=='Ro18']))
    if (num_points <= length(bra.timepoints)) {
      bra.time.overlap <- bra.timepoints[1:num_points]
    } else {
      bra.time.overlap <- bra.timepoints
    }
    atdf$sc.norm.cpm <- atdf$norm.cpm / mean(atdf$norm.cpm[atdf$timepoint %in% ara.time.overlap])        
    brdf$sc.norm.cpm <- brdf$norm.cpm / mean(brdf$norm.cpm[brdf$timepoint %in% bra.time.overlap])
    
    # get the correct times for both - move both s.t. overlap start is 0 (nb this has built in divide brassica time by 2, as 2 days between ro18 samples.)
    # if shift=5, then means the 5th bra time overlaps with the last ara time. if shift=13, then the 13th bra time overlaps the last ara time.
    # just for the sanity checking graphs! - points not actually used in the parentCSI!!

    atdf$shifted_time <- atdf$timepoint-min(ara.timepoints)
    
    br.time.key <- data.frame('time'=sort(unique(brdf$timepoint)), 'shifted_time'=0:(length(unique(bra.timepoints))-1))
    br.time.key$shifted_time <- br.time.key$shifted_time + (10-num_points)
    
    brdf <- merge(brdf, br.time.key, by.x='timepoint', by.y='time')
    
    df <- rbind(atdf, brdf)
    
    # can't have negative numbers in timepoint, as hyphens in labels!!!
    df$shifted_time <- df$shifted_time + abs(min(df$shifted_time))+1
    
    #df$shift <- num_points/3
    df$shift <- num_points
    #adf <- data.frame(accession='Col0', timepoint=at_times, sc.norm.cpm=AtVec, symbol=cur_gene, shift=num_points/3)
    #bdf <- data.frame(accession='Ro18', timepoint=br_times, sc.norm.cpm=BrVec, symbol=cur_gene, shift=num_points/3)
    #df <- rbind(adf, bdf)
    shifted_exp <- c(shifted_exp, list(df))
  }
  shifted_exp <- do.call('rbind', shifted_exp)
  
  return(shifted_exp)
}



######START OF MAIN#######
rds <- commandArgs(T)[1] # 'ro18_chiifu_apex' # 'ro18_chiifu_apex'
read_path <- commandArgs(T)[2]  # where to look for table of gene expression data
write_path <- commandArgs(T)[3]  # where to write processed data 
amalgamate <- commandArgs(T)[4] # whether should amalgamate the brassica paralogues of each arabidopsis gene.
genes_of_interest = commandArgs(T)[5:length(commandArgs(T))]

######################################################
#need to setwd() set path to correct working directory - directory with the script in)!
setwd('/jic/research-projects/bravo/alex/parent_CSI/') 
######################################################
  
  
# get Arabidopsis gene identities of interest column
ara.genes_of_interest <- paste0(tstrsplit(genes_of_interest, '-')[[1]], '-', tstrsplit(genes_of_interest, '-')[[2]])

# load table mapping arabidopsis to brassica homologues
ID_TABLE <- readRDS(paste0(read_path, "/ID_TABLE_brapa-v3.rds"))
ID_TABLE$locus_name[ID_TABLE$locus_name=='-'] <- ''



###### Load gene expression data & filter for genes of interest only
rdses <- strsplit(rds, ',')[[1]]
curr.rds <- rdses[1]
for (curr.rds in rdses) {

  path_to_rds <- paste0(read_path, curr.rds, '.rds')
  normed_expression <- readRDS(path_to_rds)
  
  # remove superfluous column
  normed_expression$est_counts <- NULL

  # filter from all genes to genes of interest only 
  if (grepl('_chiifu_', curr.rds)) { # if rds is ro18_chiifu_apex, or ro18_chiifu_leaf, or sari14_chiifu_apex/sari14_chiifu_leaf
    data <- merge(normed_expression, ID_TABLE, by.x='CDS.model', by.y='CDS.model', allow.cartesian=TRUE, all.x=T)

    #sanitise symbol for making label
    data$symbol <- gsub('/', '', data$symbol)
    data$symbol <- gsub('-', '', data$symbol)
    data$symbol <- gsub('\\.', '', data$symbol)
    # remove any duplicate rows created by the sanitising. i.e. the full data frame had two diff symbols, where '-', or '/'
    # were the only difference, and they mapped to the same gene. Causes problems as now potentially get two (identical) sets
    # of expression data for a gene of interest
    data <- unique(data)
    data$label <- paste0(data$symbol, '-', data$locus_name, '-', data$CDS.model)
    data$TAIR_description <- NULL
    # get rid of "/" to make compatible with the name sanitised list out of "parent_make_gene_list.R" - needed to make file names
    
    # cut down to only include the genes of interest - this way for amalgamating, can control the brassica copies amalgamated!
    data <- data[data$label %in% genes_of_interest, ]
    
  } else if (curr.rds == 'klepikova') {
    names(normed_expression)[names(normed_expression)=='CDS.model'] <- 'gene_id'
    data <- merge(normed_expression, ID_TABLE, by.x='gene_id', by.y='gene_id', allow.cartesian=TRUE, all.x=T) # don't want to discard any loci, just because not in ID_TABLE!
    data <- unique(data[, c('sample_id', 'accession', 'tissue', 'timepoint', 'dataset', 'group', 'norm.cpm', 'gene_id',
                     'symbol', 'locus_name')])
    data$symbol <- gsub('/', '', data$symbol)
    data$symbol <- gsub('-', '', data$symbol)
    data$symbol <- gsub('\\.', '', data$symbol)
    # remove any duplicate rows created by the sanitising. i.e. the full data frame had two diff symbols, where '-', or '/'
    # were the only difference, and they mapped to the same gene. Causes problems as now potentially get two (identical) sets
    # of expression data for a gene of interest
    data <- unique(data)
    data$locus_name <- data$gene_id # here for whatever reason, gene-Id is WAY more complete record than locus_name, but otherwise equivalent...
    data$label <- paste0(data$symbol, '-', data$locus_name) 
    # get rid of "/" to make compatible with the name sanitised list out of "parent_make_gene_list.R" - needed to make file names
    
    # cut down to only include the genes of interest - this way for amalgamating, can control the brassica copies amalgamated!
    data <- data[data$label %in% ara.genes_of_interest, ]
    
  } else {
      print('ERROR: rds is not one of ro18_chiifu_*, or klepikova!!')
      stopifnot(1==2)
    }
  
  # if want to amalgamate the expression of the paralogue genes (sum the brassicas), also need to change
  # the labels to refect the summation.
  data <- data[, c('norm.cpm', 'label', 'sample_id', 'locus_name', 'group')]
  data <- unique(data) # get rid of double counting for GO rows
  if (amalgamate=='sum') {
    print('*****************************************')
    print('Amalgamating paralogues by SUM!!')
    print('*****************************************')
    
    new.labels <- tstrsplit(data$label, '-')
    new.labels <- paste0(new.labels[[1]], '-', new.labels[[2]])
    data$new.labels <- new.labels
    data[, amalg.norm.cpm:=sum(norm.cpm), by=.(sample_id, new.labels)]
    data$label <- data$new.labels
    data$norm.cpm <- data$amalg.norm.cpm
    data$new.labels <- NULL
    data$amalg.norm.cpm <- NULL
    data <- unique(data)
  }

  label_bits <- tstrsplit(data$label, '-')
  data$ara.label <- paste0(label_bits[[1]], '-', label_bits[[2]])
  data$rds <- curr.rds
  
  if ( !exists('expression.df')) {
    expression.df <- data
  } else {
    expression.df <- rbind(expression.df, data)
  }
}

curr_expression <- expression.df

# correct for extra STM locus (AT4G37930 is Serine Transhydroxymethyltransferase 1 (STM), 
# but only want AT1G62360 (ShootMeristemless STM)!)
curr_expression <- curr_expression[curr_expression$locus_name != 'AT4G37930',] 
curr_expression <- unique(curr_expression)

rm(data, ID_TABLE, normed_expression)



##### rescale gene expression data for each gene, retaining info for each sample. 
# scale data to 0,1 for Gaussian Process step in parent_CSI.py
D <- data.table(curr_expression[, c('norm.cpm', 'label', 'sample_id', 'group', 'rds')])
D <- unique(D) # somehow some genes need this!
names(D) <- c('norm.cpm', 'gene.label', 'sample_id', 'group', 'rds')
D <- D[ , norm.cpm.sc := scale(norm.cpm) , by=.(gene.label)]

# record details for unscaling for plotting
U <- D
U$norm.cpm.sc <- NULL
U <- U[, sd:=sd(norm.cpm), by=gene.label]
U <- U[, x_bar:=mean(norm.cpm), by=gene.label]
U <- unique(U[, c('gene.label', 'x_bar', 'sd')])



##### convert to broad format 
D$norm.cpm <- NULL
D.br <- dcast(D, sample_id+group~gene.label, value.var='norm.cpm.sc')

# remove columns with any NaN - occur because genes of interest selected based on sum per ATG, so one copy of the genes
# might not be expressed - 0 to NA in rescaling
keepCols <- names(D.br)[D.br[, colSums(is.na(D.br)) != nrow(D.br)] == T]
D.br <- D.br[, keepCols, with=F]
U <- U[gene.label %in% names(D.br), ]



##### write to file 
if (!dir.exists(write_path)) {
  dir.create(write_path, recursive=TRUE)
}

fwrite(U, file=paste0(write_path, 'summary_stats.csv'))
fwrite(D.br, file=paste0(write_path, 'express_data.csv'))
