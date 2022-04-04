# script for cross validation

####################
###### SET UP ######
####################

source("00_functions.R")

# read in parameter file (specified on command line)
source("parameter_files/parameters_for_r.R")

# move to working directory 
setwd(working_dir)

## load libraries
library(moosefun)
library(tidyverse)
library(data.table)

####################################################################################
####################################################################################
## read in data
####################################################################################
####################################################################################

# model data
data_for_lm <- read.table(file=paste0(data_intermediate_dir,"06.0_data_for_lm.txt"),sep="\t",h=T)

# covar data
covar_metab <- read.table(file=paste0(data_intermediate_dir,"06.1_covar_data_for_lm.txt"), sep="\t",h=T)

# list of associated metabs
metab_list <- read.table(paste0(data_output_dir, "08.2_linear_reg_assoc_list.txt"),h=F)

# read updated feature metadata (new IDs and corrected labels included)
orig_feature <- read.table(paste0(data_intermediate_dir,"06.2_updated_orig_feature.txt"), 
                           h=T, stringsAsFactors = F,sep="\t")
dim(orig_feature)
# add feature id for merging
orig_feature$feature_id <- rownames(orig_feature)

# Table S3 (tsir results need to be added)
tableS6_prelim <- fread(paste0(data_output_dir,"06.4_linear_analysis_results_adj_model_TableS6prelim.txt"), sep = "\t")

####################################################################################
####################################################################################
## run model
####################################################################################
####################################################################################

data_for_lm$score_group <- as.factor(data_for_lm$score_group)

# set seed to ensure CV results reproducible
set.seed(63521)

# linear analysis function to include implementation of TSIR (two-step iterative resampling)
linear_analysis_tsir <- function(data, var_list, covar_metab) {
  dtst <- merge(data,covar_metab,by.x = "CLIENT_IDENTIFIER", by.y = "aln_qlet")   ## define dataset 
  linear_results = as.data.frame(matrix(data = NA, nrow = 100, ncol = nrow(var_list)))
  if (levels(dtst$score_group)[1] == "high") {
    dtst <- within(dtst,score_group <- relevel(score_group, ref = "low"))
  }
  s = base::seq(from=1,to=nrow(var_list), by=10)
  # sampling loop
  for (j in 1:100) {
    # randomly select rows for discovery
    rows_to_select <- seq(1,nrow(dtst),1)
    rows_for_disc <- sample(rows_to_select, size=525, replace = FALSE, prob = NULL)
    disc_dtst <- dtst[rows_for_disc,]
    rep_dtst <- dtst[-rows_for_disc,]
    # modelling loop
    for (i in 1:nrow(var_list)){
      if(i %in% s){
        print(paste0( "Now processing metabolite ", i, " of ", nrow(var_list) ))
      }
      mtb_name <- var_list[i,1]
      mtb_col_ref <- which(colnames(disc_dtst) == mtb_name)
      names(linear_results)[i] <- mtb_name
      # subset data to single feature
      disc_mtb <- disc_dtst[!is.na(disc_dtst[mtb_col_ref]),c("SAMPLE_NAME","score_group",mtb_name)]
      rep_mtb <- rep_dtst[!is.na(disc_dtst[mtb_col_ref]),c("SAMPLE_NAME","score_group",mtb_name)]
      # run discovery
      # check numbers (important if using unimputed data)
      group_fac <- length(unique(disc_mtb$score_group))
      if ( (group_fac > 1) & (nrow(disc_mtb) > 10) ) {
        # rank transform
        disc_mtb[,3] <- rntransform(disc_mtb[,3])
        fit <- lm(as.formula(paste0(mtb_name,"~ score_group")), 
                  data=na.omit(disc_mtb))
        coef = summary(fit)$coefficients
        group_coef = coef[2 ,]
        disc_p <- group_coef[4]
        #linear_results[j,i] <- group_coef[1]
        
      } else {
        linear_results[j,i] <- NA
      }
      # run replication
      # check numbers (important if using unimputed data)
      group_fac <- length(unique(rep_mtb$score_group))
      if ( (group_fac > 1) & (nrow(rep_mtb) > 10) ) {
        # rank transform
        rep_mtb[,3] <- rntransform(rep_mtb[,3])
        fit <- lm(as.formula(paste0(mtb_name,"~ score_group")), 
                  data=na.omit(rep_mtb))
        coef = summary(fit)$coefficients
        group_coef = coef[2 ,]
        rep_p <- group_coef[4]
        #linear_results[j,i] <- group_coef[1]
        
      } else {
        linear_results[j,i] <- NA
      }
      # input results (1=disc and rep p meet criteria)
      if (disc_p < 0.002 & rep_p < 0.05) {
        linear_results[j,i] <- 1
      }
      else {
        linear_results[j,i] <- 0
      }
    }
  }
    # save out result
    return(linear_results)
}

# run linear regression TSIR function for 29 associated metabolites
run_tsir_1 <- linear_analysis_tsir(data_for_lm, metab_list,covar_metab )

# extract per metabolite results (number of times both discovery and replication criteria met)
tsir_result <- as.data.frame(colSums(run_tsir_1))
tsir_result$feature_id <- row.names(tsir_result)
names(tsir_result)[1] <- "tsir_n_success"

# merge with metab labels
tsir_result_labelled <- left_join(tsir_result,orig_feature,by=c("feature_id"))

# add to table S5
tableS6_prelim$TSIR.result <- tsir_result_labelled$tsir_n_success[match(tableS6_prelim$Metabolite.name,tsir_result_labelled$BIOCHEMICAL)]

# write out
fwrite(tableS6_prelim, file = paste0(data_output_dir,"ForPaper/TableS6_09.1_linear_analysis_results_adj_model.txt"), sep = "\t")


#################
###### END ######
#################

print(sessionInfo())

q()