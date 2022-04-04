# script to logistic regression analysis to identify metabolite set associated with BMI score group


####################
###### SET UP ######
####################

# read in parameter file (specified on command line)
source("parameter_files/parameters_for_r.R")
#source("00_functions.R")

# move to working directory 
setwd(working_dir)

## load libraries
library(car)
library(heatmap3)
library(RColorBrewer)
library(tidyverse)
library(lme4)
library(patchwork)
library(data.table)
library(rstatix)
library(ggpubr)
library(rsq)
library(ggplot2)

####################################################################################
####################################################################################
## read in data
####################################################################################
####################################################################################

## read phenotype data from 05 script
phenofile <- read.table(file = paste0(data_intermediate_dir,"phenotype_dat_clean_with_age.txt"), 
                        h=T, stringsAsFactors = F)
phenofile$aln_qlet <- paste(phenofile$aln,phenofile$qlet,sep="_")

## read updated feature metadata (new IDs and corrected labels included)
orig_feature <- read.table(paste0(data_intermediate_dir,"06.2_updated_orig_feature.txt"), 
                           h=T, stringsAsFactors = F,sep="\t")
orig_feature$feature_id <- rownames(orig_feature)

## read xenobiotics data
load(paste0(data_intermediate_dir,"high_missingness_metab_data.RData"))
# metab_hmiss contains metabolites data
# feature_list_high_miss contains metabolites IDs

## create data frame for analysis
data <- merge(phenofile[,c("aln_qlet","pat_social_class","mat_social_class","sex","age_wks_F24","bmi_F24")], 
              metab_hmiss, by.x = "aln_qlet", by.y = "CLIENT_IDENTIFIER", all=FALSE)

###################################################
## perform logistic regression with simple model ##
###################################################

logistic_results <- as.data.frame(matrix(data = NA, nrow = length(feature_list_high_miss), ncol = ))
for (i in 1:length(feature_list_high_miss)) {
  mtb_name <- feature_list_high_miss[i]
  df_tmp <- na.omit(data[,c("aln_qlet","score_group",mtb_name)])
  if (levels(df_tmp$score_group)[1] == "high") {
    df_tmp <- within(df_tmp,score_group <- relevel(score_group, ref = "low"))
  }
  logistic_results[i,1] <- mtb_name
  logistic_results[i,2] <- nrow(df_tmp)
  logistic_results[i,3] <- table(df_tmp[,3])["1"]
  if (logistic_results[i,3] > 10) {
    logistic_results[i,4] <- "logistic"
    fit <- glm(as.formula(paste0(mtb_name,"~ score_group")),
               data = na.omit(df_tmp), family = "binomial")
    coef = summary(fit)$coefficients
    group_coef = coef[2 ,]
    eta = rsq(fit)
    logistic_results[i,5:8] <- group_coef
    logistic_results[i,9] <- group_coef[1] - (group_coef[2]*1.96)
    logistic_results[i,10] <- group_coef[1] + (group_coef[2]*1.96)
    logistic_results[i,11] <- eta
  } else {
    logistic_results[i,4:11] <- NA
  }
}
names(logistic_results)[1] <- "feature_id"
names(logistic_results)[2] <- "n_samples"
names(logistic_results)[3] <- "n_presence"
names(logistic_results)[4] <- "model_type"
names(logistic_results)[5] <- "glm_group_beta"
names(logistic_results)[6] <- "glm_group_SE"
names(logistic_results)[7] <- "glm_group_tval"
names(logistic_results)[8] <- "glm_group_p"
names(logistic_results)[9] <- "glm_group_lci"
names(logistic_results)[10] <- "glm_group_uci"
names(logistic_results)[11] <- "glm_group_r2"

## adjust p-val
logistic_results$glm_group_adjp = p.adjust(logistic_results[, "glm_group_p"], method = "BH")

## identify significant features
logistic_results$logistic_flag = 0
logistic_results[which(logistic_results$glm_group_adjp <= 0.05),"logistic_flag"] <- 1
sum(logistic_results$logistic_flag)

# sort by p
logistic_results <- logistic_results[order(logistic_results$glm_group_p,decreasing = F),]

# add rank
logistic_results$logistic_score_p_rank <- rank(logistic_results$glm_group_p)


glm_simple_res <- merge(logistic_results, orig_feature, by = "feature_id", all.x = T)
glm_simple_res <- glm_simple_res[order(glm_simple_res$logistic_score_p_rank,decreasing = F),]

# write out in full
fwrite(glm_simple_res, paste0(data_output_dir, "07.1_logistic_regression_results_simple_model.txt"), sep = "\t")

# write out results for supp table
# Paper Ref: Table S5

# tidy up column headers
names(glm_simple_res)[c(16:18,21,22:27,2,3,5:6,9,10,8,11:14)] <-
  c("Metabolite.name","Super.pathway","Sub.pathway","Metabolon.chem.id","Retention.index","Metabolite.mass","Pubchem.id","CAS.number","Kegg.id","HMDB.id",
      "Sample.size","Number.of.samples.with.measurement","Beta","Beta.standard.error","Lower.95%.CI","Upper.95%.CI","P.value",
      "VE.by.score.group","Benjamini.Hochberg.p.value","Associated.flag","Associated.rank")

fwrite(glm_simple_res[!is.na(glm_simple_res$model_type),c(16:18,21,22:27,2,3,5:6,9,10,8,11:14)], paste0(data_output_dir, "ForPaper/TableS5_07.1_logistic_regression_results_simple_model.txt"), sep = "\t")


####################################################
## perform logistic regression with complex model ##
####################################################

logistic_results <- as.data.frame(matrix(data = NA, nrow = length(feature_list_high_miss), ncol = ))
for (i in 1:length(feature_list_high_miss)) {
  mtb_name <- feature_list_high_miss[i]
  df_tmp <- na.omit(data[,c("aln_qlet","score_group",mtb_name,"pat_social_class","mat_social_class")])
  if (levels(df_tmp$score_group)[1] == "high") {
    df_tmp <- within(df_tmp,score_group <- relevel(score_group, ref = "low"))
  }
  logistic_results[i,1] <- mtb_name
  logistic_results[i,2] <- nrow(df_tmp)
  logistic_results[i,3] <- ifelse(table(df_tmp[,3])["0"] != nrow(df_tmp), table(df_tmp[,3])["1"], 0)
  
  if (logistic_results[i,3] > 10) {
    logistic_results[i,4] <- "logistic"
    fit <- glm(as.formula(paste0(mtb_name,"~ score_group + pat_social_class + mat_social_class")),
               data = na.omit(df_tmp), family = "binomial")
    coef = summary(fit)$coefficients
    group_coef = coef[2 ,]
    eta = rsq.partial(fit)
    logistic_results[i,5:8] <- group_coef
    logistic_results[i,9:11] <- eta$partial.rsq
  } else {
    logistic_results[i,4:11] <- NA
  }
}
names(logistic_results)[1] <- "feature_id"
names(logistic_results)[2] <- "n_samples"
names(logistic_results)[3] <- "n_presence"
names(logistic_results)[4] <- "model_type"
names(logistic_results)[5] <- "glm_group_beta"
names(logistic_results)[6] <- "glm_group_SE"
names(logistic_results)[7] <- "glm_group_tval"
names(logistic_results)[8] <- "glm_group_p"
names(logistic_results)[9:11] = paste0("r2_", eta$variable)

## adjust p-val
logistic_results$glm_group_adjp = p.adjust(logistic_results[, "glm_group_p"], method = "BH")

## identify significant features
logistic_results$logistic_flag = 0
logistic_results[which(logistic_results$glm_group_adjp <= 0.05),"logistic_flag"] <- 1
sum(logistic_results$logistic_flag)

# sort by p
logistic_results <- logistic_results[order(logistic_results$glm_group_p,decreasing = F),]

# add rank
logistic_results$logistic_score_p_rank <- rank(logistic_results$glm_group_p)

glm_complex_res <- merge(logistic_results, orig_feature, by = "feature_id", all.x = T)
glm_complex_res <- glm_complex_res[order(glm_complex_res$logistic_score_p_rank,decreasing = F),]

fwrite(glm_complex_res, paste0(data_output_dir, "07.2_logistic_regression_results_complex_model.txt"), sep = "\t")

test <- merge(glm_simple_res, glm_complex_res, by = "feature_id")
test <- test[order(test$Associated.rank),]

# write out in full
fwrite(test, paste0(data_output_dir, "07.3_logistic_regression_results_models_comparison.txt"), sep = "\t")

#################
###### END ######
#################

print(sessionInfo())

q()