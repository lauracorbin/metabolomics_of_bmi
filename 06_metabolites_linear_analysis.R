# script for linear regression analysis to identify metabolite set associated with BMI score group

####################
###### SET UP ######
####################

# read in function file
source("00_functions.R")

# read in parameter file (specified on command line)
source("parameter_files/parameters_for_r.R")

# move to working directory 
setwd(working_dir)

## load libraries
library(car)
library(heatmap3)
library(RColorBrewer)
library(plyr)
library(tidyverse)
library(lme4)
library(patchwork)
library(data.table)
library(rstatix)
library(ggpubr)
library(missForest)
library(moosefun)
library(ggplot2)
library(ggpubr)
library(ggrepel)


####################################################################################
####################################################################################
## read in data
####################################################################################
####################################################################################

## Metabolon data ##
## original scale version (after QC)
orig_data <- read.table(paste0(metab_input_dir,"qc_data/B3194_BMI_RbG_2020_02_03_QCd_data.txt"), h=T)
dim(orig_data)
head(names(orig_data))
tail(names(orig_data))

## feature metadata
orig_feature <- read.table(paste0(metab_input_dir,"B3194_BMI_RbG_2020_02_03_Metabolon_featuredata.txt"), 
                           h=T, stringsAsFactors = F)
dim(orig_feature)

## sample metadata with group added
orig_sample <- read.table(paste0(data_intermediate_dir,"Metabolon_sampledata_with_group.txt"), 
                          h=T, stringsAsFactors = F)
dim(orig_sample)

## pheno data ##
## extract from 05 script
phenofile <- read.table(file = paste0(data_intermediate_dir,"phenotype_dat_clean_with_age.txt"), 
                        h=T, stringsAsFactors = F)          

## read in list of ALSPAC IDs to remove 
individuals_to_remove <- readRDS(file = paste0(data_intermediate_dir, "individuals_to_remove.rds"))

## read in extra chemical IDs received from Metabolon. I sent a list of X IDs I was interested in (ie those that appeared in the 'hits' list)
## Metabolon sent these back in Feb 2020
extra_biochem_ids <- read.table(file = paste0(data_intermediate_dir,"extra_biochem_ids.txt"), 
                                  h=T, stringsAsFactors = F, sep="\t") 

# update biochem IDs to include new info for associated
for (i in 1:nrow(extra_biochem_ids)) {
  # find columns to replace
  biochem_col <- which(names(orig_feature) == "BIOCHEMICAL")
  super_col <- which(names(orig_feature) == "SUPER_PATHWAY")
  sub_col <- which(names(orig_feature) == "SUB_PATHWAY")
  # look for biochem
  orig_feature_row <- which(orig_feature$BIOCHEMICAL == extra_biochem_ids[i,1])
  # replace ids
  orig_feature[orig_feature_row,biochem_col] = extra_biochem_ids[i,2]
  orig_feature[orig_feature_row,super_col] = extra_biochem_ids[i,3]
  orig_feature[orig_feature_row,sub_col] = extra_biochem_ids[i,4]
}

## update biochem IDs to include corrected labels

# read in information relating to a metabolite label correction issued by Metabolon in Jan 2022.
new_labels <- read.table(file="input/metabolon_id_update.csv",sep=",",h=T)

for (i in 1:nrow(new_labels)) {
  # find column
  to_replace <- which(orig_feature$BIOCHEMICAL == new_labels[i,c("IncorrectID")])
  if (length(to_replace) > 0){
    print(paste0("replacing ", new_labels[i,2], " with ", new_labels[i,3]))
    # chem id
    orig_feature[to_replace,c("CHEMICAL_ID")] <- new_labels[i,4]
    # comp id - don't update as is used to merge with metab data (not included in output anyway)
    #orig_feature[to_replace,c("COMP_ID")] <- new_labels[i,5]
    # CAS
    orig_feature[to_replace,c("CAS")] <- new_labels[i,6]
    # pubchem
    orig_feature[to_replace,c("PUBCHEM")] <- new_labels[i,7]
    # HMDB
    orig_feature[to_replace,c("HMDB")] <- new_labels[i,8]
    # name
    orig_feature[to_replace,c("BIOCHEMICAL")] <- new_labels[i,3]
    # set missing where no info
    orig_feature[to_replace,c("KEGG")] <- NA
    # update row name - don't update as is used to merge with metab data (not included in output anyway)
    #row.names(orig_feature)[to_replace] <- paste0("compid_",orig_feature[to_replace,"COMP_ID"])
  }
  else {
    print(paste0("Skipping ", new_labels[i,c("IncorrectID")]))
  }
}
dim(orig_feature)

## save out updated feature info
write.table(orig_feature,
            file=paste0(data_intermediate_dir,"06.2_updated_orig_feature.txt"),
            quote=T, sep="\t")

####################################################################################
####################################################################################
## pheno process and checks
####################################################################################
####################################################################################

# check group allocation
table(orig_sample$score_group)

# add aln_qlet ID to phenofile (for use with merging)
phenofile$aln_qlet <- paste(phenofile$aln,phenofile$qlet,sep="_")
names(phenofile)

phenofile <- subset(phenofile, select = -c(score_group))

####################################################################################
####################################################################################
## make data files for running models on
####################################################################################
####################################################################################

## generate feature list
feature_list <- names(orig_data)[1:(ncol(orig_data))]
head(feature_list)
length(feature_list)

# bring in sample id
orig_data$SAMPLE_NAME <- rownames(orig_data)

# drop duplicates from metabolite file based on previously  determined criteria (sample missingness)
sample_to_remove <- readRDS(file = paste0(data_intermediate_dir, "samples_to_remove.rds"))
metab_nodups <- orig_sample[-which(orig_sample$SAMPLE_NAME %in% sample_to_remove),]
dim(metab_nodups)
table(metab_nodups$score_group)

# merge, keeping only those that are in both files
metab_full <- merge(metab_nodups,orig_data,by="SAMPLE_NAME", all=FALSE)
dim(metab_full)

# change group to factor
metab_full$score_group <- as.factor(metab_full$score_group)
rownames(metab_full) <- metab_full$CLIENT_IDENTIFIER

# drop those to remove
dim(metab_full)
metab_full <- metab_full[-which(metab_full$CLIENT_IDENTIFIER %in% individuals_to_remove),]
dim(metab_full)
table(metab_full$score_group)

####################################################################################
####################################################################################
## Tests for batch effects
####################################################################################
####################################################################################

dtst <- metab_full

# fisher test: BOX_ID/BOX_WELL/RUN_DAY ~ score_group
dtst$sample_id <- rownames(dtst)
dat_box_res <- as.data.frame(matrix(data = NA, nrow = 3, ncol = 3))
for (i in 1:length(c("BOX_ID","BOX_WELL","RUN_DAY"))) {
  var_name <- c("BOX_ID","BOX_WELL","RUN_DAY")[i]
  df_tmp <- na.omit(dtst[,c(var_name,"score_group")])
  df_tmp[,var_name] <- factor(df_tmp[,var_name])
  res <- fisher.test(table(df_tmp[[var_name]],df_tmp$score_group),simulate.p.value = TRUE)
  dat_box_res[i,1:3] <- c(var_name, res$p.value, dim(df_tmp)[1])
}
names(dat_box_res) <- c("variable", "pval", "N")
write.table(dat_box_res,
            file=paste0(data_output_dir,"06.1_metabolites_batch_effect_analysis_raw.txt"),
            quote=F, sep="\t",row.names=F)


####################################################################################
####################################################################################
## QC by data missingness
####################################################################################
####################################################################################

# the percentage of missing data of metabolites
metab_miss_percent <- apply(metab_full[,feature_list],2,function(x) sum(is.na(x))/length(x))
sample_miss_percent <- apply(metab_full[,feature_list],1,function(x) sum(is.na(x))/length(x))
summary(metab_miss_percent)
summary(sample_miss_percent)

p_metab_miss <- ggplot() +
  geom_histogram(aes(metab_miss_percent),bins = 30) +
  labs(title = "features") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))

p_sample_miss <- ggplot() +
  geom_histogram(aes(sample_miss_percent),bins = 30) +
  labs(title = "samples") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))

p_metab_miss | p_sample_miss

ggsave(path = data_output_dir,filename = "06.2_metab_data_missingness.png",
       plot = last_plot(),width = 20,height = 10,units = c("cm"),dpi = 320)

# group metabolites based on missingness
missingness_threshold <- 0.2
feature_list_high_miss <- feature_list[feature_list %in% names(which(colMeans(is.na(metab_full[,feature_list])) > missingness_threshold))]
print(paste0("There are ", length(feature_list_high_miss), " metabolites with more than 20% missing."))
feature_list_low_miss <- feature_list[!feature_list %in% names(which(colMeans(is.na(metab_full[,feature_list])) > missingness_threshold))]
print(paste0("There are ", length(feature_list_low_miss), " metabolites with less than or equal to 20% missing."))

# calculate metabolites missingness by group
missingness_by_group <- metab_full[c(feature_list_low_miss,"score_group")] %>% 
  group_by(score_group) %>%
  dplyr::summarise(across(starts_with("compid_"), .fns = is.na)) %>%
  dplyr::summarise(across(starts_with("compid_"), .fns = mean)) %>%
  gather(variables, val, -score_group) %>%
  spread(score_group, val)

# convert metabolites with high missingness to presence/absence
# because previous qc has removed non-xenobiotics metabolites with missingness greater than 20%
metab_full[,feature_list_high_miss] <- pa.convert(metab_full[,feature_list_high_miss])
metab_hmiss <- metab_full[,c(names(metab_full)[1:11],feature_list_high_miss)]
save(metab_hmiss,feature_list_high_miss, 
     file = paste0(data_intermediate_dir, "high_missingness_metab_data.RData"))

## impute missing data for metabolites with low missingness
dat_to_impute <- metab_full[,feature_list_low_miss]

# random forest
set.seed(63525)
#mf <- missForest(dat_to_impute)     # this takes about 1 hour to complete
#save(mf,file = paste0(data_intermediate_dir, "imputation_missForest_result.RData"))
load(paste0(data_intermediate_dir, "imputation_missForest_result.RData"))
imputed_dat_mf <- as.data.frame(mf$ximp)

# replace missing data with imputed values
metab_imputed_mf <- metab_full
metab_imputed_mf[,feature_list_low_miss] <- imputed_dat_mf

# data summary for raw dataset
dtst <- metab_full[,feature_list]
ds_results <- as.data.frame(t(apply(dtst,2,function(x) { 
  o = summary(x)
  if(length(o) == 6){
    o = c(o, 0)
    names(o)[7] = "NA's"
  }
  return(o)
})))
write.table(ds_results,
            file=paste0(data_output_dir,"06.3_metabolite_summaries_raw.txt"),
            quote=F, sep="\t",row.names=F)

# data summary for random forest imputed dataset
dtst <- metab_imputed_mf[,feature_list]
ds_results <- as.data.frame(t(apply(dtst,2,function(x) { 
  o = summary(x)
  if(length(o) == 6){
    o = c(o, 0)
    names(o)[7] = "NA's"
  }
  return(o)
})))
write.table(ds_results,
            file=paste0(data_output_dir,"06.3_metabolite_summaries_imputed.txt"),
            quote=F, sep="\t",row.names=F)

####################################################################################
####################################################################################
## run linear model (low missingness metabolites)
####################################################################################
####################################################################################

data_for_lm <- metab_imputed_mf

# save out for use in cross validation script
write.table(data_for_lm,
            file=paste0(data_intermediate_dir,"06.0_data_for_lm.txt"),
            quote=F, sep="\t",row.names=F)

## define covariates and other data
var_to_extract <- c("aln_qlet","mat_social_class","pat_social_class","bmi_F24","sex","age_wks_F24")
covar_metab <- phenofile[,grepl(paste(var_to_extract, collapse="|"), names(phenofile))]
covar_metab$mat_social_class <- factor(covar_metab$mat_social_class)
covar_metab$pat_social_class <- factor(covar_metab$pat_social_class)

# save out for use in cross validation script
write.table(covar_metab,
            file=paste0(data_intermediate_dir,"06.1_covar_data_for_lm.txt"),
            quote=F, sep="\t",row.names=F)

res_simple <- linear_analysis_simple(data = data_for_lm, var_list = feature_list_low_miss, covar_metab = covar_metab)
res_complex <- linear_analysis_complex(data = data_for_lm, var_list = feature_list_low_miss, covar_metab = covar_metab) # sensitivity analysis
res_BMI <- linear_analysis_BMI(data = data_for_lm, var_list = feature_list_low_miss, covar_metab = covar_metab) # sensitivity analysis

# Q-Q plot of p-values
plots_simple <- linear_analysis_plots(data = res_simple, model = "simple")
ggsave(path = data_output_dir,
       filename = paste0("06.5_histogram_rawpval_metabolite_simple_model.png"),
       plot = plots_simple[[1]], width = 15, height = 15, units = c("cm"), dpi = 320)
ggsave(path = data_output_dir,
       filename = paste0("06.5_qq_plot_rawpval_metabolite_simple_model.png"),
       plot = plots_simple[[2]], width = 15, height = 15, units = c("cm"), dpi = 320)

plots_complex <- linear_analysis_plots(data = res_complex, model = "complex")
ggsave(path = data_output_dir,
       filename = paste0("06.5_histogram_rawpval_metabolite_complex_model.png"),
       plot = plots_complex[[1]], width = 15, height = 15, units = c("cm"), dpi = 320)
ggsave(path = data_output_dir,
       filename = paste0("06.5_qq_plot_rawpval_metabolite_complex_model.png"),
       plot = plots_complex[[2]], width = 15, height = 15, units = c("cm"), dpi = 320)

plots_BMI <- linear_analysis_plots_BMI(data = res_BMI, model = "BMI")
ggsave(path = data_output_dir,
       filename = paste0("06.5_histogram_rawpval_metabolite_BMI_model.png"),
       plot = plots_BMI[[1]], width = 15, height = 15, units = c("cm"), dpi = 320)
ggsave(path = data_output_dir,
       filename = paste0("06.5_qq_plot_rawpval_metabolite_BMI_model.png"),
       plot = plots_BMI[[2]], width = 15, height = 15, units = c("cm"), dpi = 320)

# combine metabolites information with linear results
orig_feature$feature_id <- rownames(orig_feature)
merged_res <- merge(res_simple, orig_feature, by = "feature_id", all.x = T)

# write out in full
fwrite(merged_res, file = paste0(data_output_dir,"06.4_linear_analysis_results_simple_model.txt"), sep = "\t")

# write out results for supp table
# add missingness info
metab_miss_percent_table <- as.data.frame(metab_miss_percent)
metab_miss_percent_table$feature_id <- rownames(metab_miss_percent_table)
merged_res$prop_missing <- metab_miss_percent_table$metab_miss_percent[match(merged_res$feature_id,metab_miss_percent_table$feature_id)]
merged_res$percent_missing <- (merged_res$prop_missing)*100 

## calculate fold change on raw data to add ot output table
fold.change.out <- as.data.frame(matrix(data=NA,nrow=(ncol(metab_full)-11),ncol=4))
for (i in 12:ncol(metab_full)) {
  out_row <- i-11
  metab_name <- names(metab_full[i])
  # subset data to single feature
  mtb <- metab_full[!is.na(metab_full[i]), c("score_group",metab_name)]
  # calculate means by group
  median_low <- median(mtb[mtb$score_group == "low",2])
  median_high <- median(mtb[mtb$score_group == "high",2])
  # calculate fold change and log2 of FC
  fc <- median_high/median_low
  # populate output daataframe
  fold.change.out[out_row,1] <- metab_name
  fold.change.out[out_row,2] <- median_low
  fold.change.out[out_row,3] <- median_high
  fold.change.out[out_row,4] <- log2(fc)
}
names(fold.change.out)[1] <- "feature_id"
names(fold.change.out)[2] <- "median_abundance_low_bmi"
names(fold.change.out)[3] <- "median_abundance_high_bmi"
names(fold.change.out)[4] <- "log2_fc_high_over_low_bmi"

merged_res$log2_fc_high_over_low_bmi <- fold.change.out$log2_fc_high_over_low_bmi[match(merged_res$feature_id,fold.change.out$feature_id)]

# replace sub-pathway NA with unknown
merged_res <- merged_res %>%
  mutate(SUPER_PATHWAY = coalesce(SUPER_PATHWAY, "Unknown"),
         SUB_PATHWAY = coalesce(SUB_PATHWAY, "Unknown"))

########
# Volcano plot [code from David Hughes]
# Paper ref: Figure 4

## define working data
mydata = data.frame(  beta = merged_res$log2_fc_high_over_low_bmi, p = -log10(merged_res$lm_group_p) , cat = merged_res$SUPER_PATHWAY,
                        adjp = merged_res$lm_group_adjp, biochemical = merged_res$BIOCHEMICAL)

## add in labelling var
mydata$ToLabel <- 0
## labels for representative associated (list derived from code applied later in script 08)
assoc_rep <- c("metabolonic lactone sulfate","cortisone","sphingomyelin (d18:2/16:0, d18:1/16:1)*","tridecenedioate (C13:1-DC)*",
              "perfluorooctanesulfonate (PFOS)","hippurate","bilirubin degradation product, C16H18N2O5 (1)*",
              "3-hydroxy-2-ethylpropionate","pregnenolone sulfate","O-sulfo-L-tyrosine","glycocholenate sulfate*")

mydata[which(mydata$biochemical %in% assoc_rep),c("ToLabel")] <- 1

##PLOT
pdf_out = paste0(data_output_dir,"ForPaper/Fig4_06.1_volcano.pdf")
pdf(pdf_out)

## volcano plot
volcano = mydata %>% ggplot(aes(x = beta, y = p)) +
  geom_point( aes(color = cat)) + 
  scale_color_brewer(palette="Spectral") + 
  #scale_color_manual(values=brewer.pal(11,"Paired")) +
  labs(x = "log2 median fold change (high/low)", y = "-log10(P)") +
  geom_density_2d( color = "grey50", show.legend=FALSE) +
  labs(size = "Sample Size", color = paste0("Super-pathway") ) +
  guides(col = guide_legend(ncol = 3, byrow = TRUE, override.aes = list(size=1) )) +
  geom_hline(yintercept=-log10(0.0016), col="black",lty=2) +
  theme_minimal() +
  guides(shape = guide_legend(override.aes = list(size = 0.5))) +
  geom_text_repel(aes(label=ifelse(ToLabel == 1, as.character(biochemical),'')),hjust=1.2,vjust=0.5,cex=2.5) +
  theme(legend.position = "bottom",plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))
volcano

invisible(dev.off())

# Paper Ref: Table S4

# tidy up column headers
names(merged_res)[c(17:19,22,23:28,30,5:6,9:10,8,4,13:15,31)] <-
  c("Metabolite.name","Super.pathway","Sub.pathway","Metabolon.chem.id","Retention.index","Metabolite.mass","Pubchem.id","CAS.number","Kegg.id","HMDB.id",
      "%.missing.values","Beta","Beta.standard.error","Lower.95%.CI","Upper.95%.CI","P.value","Model.R2","Benjamini.Hochberg.p.value",
      "Associated.flag","Associated.rank","Log2.fold.change.(high.cf.low)")
fwrite(merged_res[,c(17:19,22,23:28,30,5:6,9:10,8,4,13:15,31)], file = paste0(data_output_dir,"ForPaper/TableS4_06.4_linear_analysis_results_simple_model.txt"), sep = "\t")

########
# combine metabolites information with linear results (fully adjusted model)
merged_res_complex <- merge(res_complex, orig_feature, by = "feature_id", all.x = T)
merged_res_complex <- merged_res_complex[order(merged_res_complex$lm_group_p),]

# write out results for supp table
# Paper Ref: Table S6 prelim (info added in script 09)

## ID those metab assoc in main 
associated_list <- merged_res[merged_res$Associated.flag == 1, "Metabolite.name"]

# tidy up column headers
names(merged_res_complex)[c(19,5:6,9,10,8,11:13,4,15)] <- 
  c("Metabolite.name","Beta","Beta.standard.error","Lower.95%.CI","Upper.95%.CI","P.value","VE.by.score.group",
      "VE.by.paternal.social.class","VE.by.maternal.social.class","Model.R2","Benjamini.Hochberg.p.value")

fwrite(merged_res_complex[merged_res_complex$Metabolite.name %in% associated_list,c(19,5:6,9,10,8,11:13,4,15)], file = paste0(data_output_dir,"06.4_linear_analysis_results_adj_model_TableS6prelim.txt"), sep = "\t")

########
# combine metabolites information with linear results (bmi analysis)
merged_res_bmi <- merge(res_BMI, orig_feature, by = "feature_id", all.x = T)
merged_res_bmi <- merged_res_bmi[order(merged_res_bmi$lm_BMI_p),]

########
# combine three sets of results to write out
names(res_complex)[2:(dim(res_complex)[2])] <- paste0("c_",names(res_complex)[2:(dim(res_complex)[2])])
names(res_BMI)[2:(dim(res_BMI)[2])] <- paste0("BMI_",names(res_BMI)[2:(dim(res_BMI)[2])])
merged_res <- merge(merged_res, res_complex, by = "feature_id")
merged_res <- merge(merged_res, res_BMI, by = "feature_id")

merged_res <- merged_res[order(merged_res$P.value),]

# write out in full
fwrite(merged_res, file = paste0(data_output_dir,"06.4_linear_analysis_results_all_models.txt"), sep = "\t")

## check correlation between primary result and that with extra adjustments
cor.test(merged_res$Beta,merged_res$c_lm_group_beta)
# restrict to associated
merge_res_assoc <- merged_res[merged_res$Associated.flag == 1,]

cor.test(merge_res_assoc$Beta,merge_res_assoc$c_lm_group_beta)

#################
###### END ######
#################

print(sessionInfo())

q()
  