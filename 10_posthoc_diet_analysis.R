# script for post hoc analysis of dietary pheno

####################
###### SET UP ######
####################

# read in parameter file (specified on command line)
source("parameter_files/parameters_for_r.R")

# move to working directory 
setwd(working_dir)

## load libraries
library(tidyverse)
library(data.table)
library(dplyr)
####################################################################################
####################################################################################
## read in data
####################################################################################
####################################################################################

## pheno data ##
## extract from 05 script
phenofile <- read.table(file = paste0(data_intermediate_dir,"phenotype_dat_clean_with_age.txt"), 
                        h=T, stringsAsFactors = F)          

## original scale version (after QC)
orig_data <- read.table(paste0(metab_input_dir,"qc_data/B3194_BMI_RbG_2020_02_03_QCd_data.txt"), h=T)

## read updated feature metadata (new IDs and corrected labels included)
orig_feature <- read.table(paste0(data_intermediate_dir,"06.2_updated_orig_feature.txt"), 
                           h=T, stringsAsFactors = F,sep="\t")
orig_feature$feature_id <- rownames(orig_feature)

## sample metadata with group added
orig_sample <- read.table(paste0(data_intermediate_dir,"Metabolon_sampledata_with_group.txt"), 
                          h=T, stringsAsFactors = F)
dim(orig_sample)

####################################################################################
####################################################################################
## derive new diet variables
####################################################################################
####################################################################################

### food pref from YPE questionnaire (age 25 years)

# fish (sum of: fried/battered, baked/steamed, prawns, salmon, shellfish, smoked salmon, tuna)
phenofile$fish.preference <- rowSums(phenofile[c("YPE9025","YPE9026","YPE9027","YPE9028","YPE9029","YPE9030","YPE9031")], na.rm=FALSE)
summary(phenofile$fish.preference)

# fresh veg (sum of: garlic, green olives, mushrooms, onions, tomoatoes, chilli peppers, artichokes, asparagus, aubergines, 
# avocados, black olives, broad beans, broccoli, brussels sprouts, cabbage, carrots, spinach)
phenofile$vegetable.preference <- rowSums(phenofile[c("YPE9012","YPE9013","YPE9014","YPE9015","YPE9075","YPE9076","YPE9086","YPE9087","YPE9088",
                                           "YPE9089","YPE9090","YPE9091","YPE9092","YPE9093","YPE9094","YPE9095","YPE9096")], na.rm=FALSE)
summary(phenofile$vegetable.preference)

# fresh fruit (sum of: apples, bananas, cherries, dried fruit, lemons, oranges, pears, strawberries)
phenofile$fruit.preference <- rowSums(phenofile[c("YPE9032","YPE9033","YPE9034","YPE9035","YPE9036","YPE9037","YPE9038","YPE9039")], na.rm=FALSE)
summary(phenofile$fruit.preference)

####################################################################################
####################################################################################
## summarise by recall group
####################################################################################
####################################################################################

# variables to summarise
char_var_list <- c("fish.preference","vegetable.preference","fruit.preference")

# make results output
char_var_summ <- as.data.frame( matrix(data = NA, nrow = length(char_var_list), ncol = 22) )

# run analysis loop
for (i in 1:length(char_var_list)) {
  print(paste0("Processing summary for: ", char_var_list[i]))
  # var name
  char_var_summ[i,1] <- char_var_list[i]

  # RbG cohort (phenofile)
  # var col
  col_num <- which(names(phenofile) == char_var_list[i])
  # RbG cohort - all
  # calculate N 
  char_var_summ[i,2] <- length(which(!is.na(phenofile[,col_num])))
  # calculate means
  char_var_summ[i,3] <- mean(phenofile[,col_num], na.rm=T)
  # calculate SD
  char_var_summ[i,4] <- sd(phenofile[,col_num], na.rm=T)
  # RbG cohort - males - high
  # calculate N 
  char_var_summ[i,5] <- length(which(!is.na(phenofile[phenofile$sex == 1 & phenofile$score_group == "high",col_num])))
  # calculate means
  char_var_summ[i,6] <- mean(phenofile[phenofile$sex == 1 & phenofile$score_group == "high",col_num], na.rm=T)
  # calculate SD
  char_var_summ[i,7] <- sd(phenofile[phenofile$sex == 1 & phenofile$score_group == "high",col_num], na.rm=T)
  # RbG cohort - males - low
  # calculate N 
  char_var_summ[i,8] <- length(which(!is.na(phenofile[phenofile$sex == 1 & phenofile$score_group == "low",col_num])))
  # calculate means
  char_var_summ[i,9] <- mean(phenofile[phenofile$sex == 1 & phenofile$score_group == "low",col_num], na.rm=T)
  # calculate SD
  char_var_summ[i,10] <- sd(phenofile[phenofile$sex == 1 & phenofile$score_group == "low",col_num], na.rm=T)
  # RbG cohort - females - high
  # calculate N 
  char_var_summ[i,11] <- length(which(!is.na(phenofile[phenofile$sex == 2 & phenofile$score_group == "high",col_num])))
  # calculate means
  char_var_summ[i,12] <- mean(phenofile[phenofile$sex == 2 & phenofile$score_group == "high",col_num], na.rm=T)
  # calculate SD
  char_var_summ[i,13] <- sd(phenofile[phenofile$sex == 2 & phenofile$score_group == "high",col_num], na.rm=T)
  # RbG cohort - females - low
  # calculate N 
  char_var_summ[i,14] <- length(which(!is.na(phenofile[phenofile$sex == 2 & phenofile$score_group == "low",col_num])))
  # calculate means
  char_var_summ[i,15] <- mean(phenofile[phenofile$sex == 2 & phenofile$score_group == "low",col_num], na.rm=T)
  # calculate SD
  char_var_summ[i,16] <- sd(phenofile[phenofile$sex == 2 & phenofile$score_group == "low",col_num], na.rm=T)
  # RbG cohort - high
  # calculate N 
  char_var_summ[i,17] <- length(which(!is.na(phenofile[phenofile$score_group == "high",col_num])))
  # calculate means
  char_var_summ[i,18] <- mean(phenofile[phenofile$score_group == "high",col_num], na.rm=T)
  # calculate SD
  char_var_summ[i,19] <- sd(phenofile[phenofile$score_group == "high",col_num], na.rm=T)
  # RbG cohort - low
  # calculate N 
  char_var_summ[i,20] <- length(which(!is.na(phenofile[phenofile$score_group == "low",col_num])))
  # calculate means
  char_var_summ[i,21] <- mean(phenofile[phenofile$score_group == "low",col_num], na.rm=T)
  # calculate SD
  char_var_summ[i,22] <- sd(phenofile[phenofile$score_group == "low",col_num], na.rm=T)
  
}

names(char_var_summ) <- c("variable","rbg_all_n","rbg_all_mean","rbg_all_sd",
                          "rbg_males_high_n","rbg_males_high_mean","rbg_males_high_sd","rbg_males_low_n","rbg_males_low_mean","rbg_males_low_sd",
                          "rbg_females_high_n","rbg_females_high_mean","rbg_females_high_sd","rbg_females_low_n","rbg_females_low_mean","rbg_females_low_sd",
                          "rbg_high_n","rbg_high_mean","rbg_high_sd","rbg_low_n","rbg_low_mean","rbg_low_sd")


####################################################################################
####################################################################################
## test for difference by recall group
####################################################################################
####################################################################################

dtst <- phenofile[,c(char_var_list, "score_group")]

# transform score group into factor
dtst$score_group <- as.factor(dtst$score_group)
dtst$score_group <- relevel(dtst$score_group, ref = "low")

# perform rank-test on phenotype measures
res_mean_difference <- data.frame()

for (i in 1:length(char_var_list)) {
  var_name <- char_var_list[i]
  res_mean_difference[i,1] <- var_name
  df_tmp <- na.omit(dtst[c(var_name, "score_group")])
  if (levels(df_tmp$score_group)[1] == "low") {
    df_tmp <- within(df_tmp,score_group <- relevel(score_group, ref = "high"))
  }
  # check distribution
  res_wilcox <- wilcox.test(as.formula(paste(var_name, "score_group", sep = "~")),
                data = df_tmp)
  res_mean_difference[i,2] <- c(res_wilcox$p.value)
  res_mean_difference[i,3] <- nrow(df_tmp)
}
colnames(res_mean_difference) <- c("phenotype","wilcox_pval","N")


## join test stats to summary tables from above and write out
char_var_summ_plus <- left_join(char_var_summ,res_mean_difference,by = c("variable" = "phenotype"))

# save out
fwrite(char_var_summ_plus, paste0(data_output_dir, "ForPaper/TableS8_10.1_diet_characterisation.txt"), sep = "\t")


####################################################################################
####################################################################################
## check for associations of food pref with specific metabs
####################################################################################
####################################################################################

### add hippurate and PFOS to phenofile
# extract compid
hippurate_compid <- orig_feature[which(orig_feature$BIOCHEMICAL == "hippurate"),c("feature_id")]
pfos_compid <- orig_feature[which(orig_feature$BIOCHEMICAL == "perfluorooctanesulfonate (PFOS)"),c("feature_id")]

# drop duplicates from metabolite file based on previously  determined criteria (sample missingness)
sample_to_remove <- readRDS(file = paste0(data_intermediate_dir, "samples_to_remove.rds"))
metab_nodups <- orig_sample[-which(orig_sample$SAMPLE_NAME %in% sample_to_remove),]
dim(metab_nodups)
table(metab_nodups$score_group)

# add IDs to metab file to enable merging with phenofile
orig_data$SAMPLE_NAME <- row.names(orig_data)
orig_data$aln_qlet <- metab_nodups$CLIENT_IDENTIFIER[match(orig_data$SAMPLE_NAME,metab_nodups$SAMPLE_NAME)]

# merge metab in with pheno
phenofile_with_metab <- left_join(phenofile,orig_data,by="aln_qlet")

# transform score group into factor
phenofile_with_metab$score_group <- as.factor(phenofile_with_metab$score_group)
phenofile_with_metab$score_group <- relevel(phenofile_with_metab$score_group, ref = "low")

# hippurate and fruit and veg
hippurate_col <- which(names(phenofile_with_metab) == hippurate_compid)
names(phenofile_with_metab)[hippurate_col] <- "hippurate"
# standardise
phenofile_with_metab$hippurate_z <- scale(phenofile_with_metab$hippurate)
# log because of skewed distribution
phenofile_with_metab$hippurate_log <- log(phenofile_with_metab$hippurate)

fit_hip_veg_pref <- lm(phenofile_with_metab$hippurate_log ~ phenofile_with_metab$vegetable.preference)
summary(fit_hip_veg_pref)
fit_hip_veg_pref <- lm(phenofile_with_metab$hippurate_log ~ phenofile_with_metab$vegetable.preference + phenofile_with_metab$sex)
summary(fit_hip_veg_pref)
fit_hip_veg_pref <- lm(phenofile_with_metab$hippurate_log ~ phenofile_with_metab$vegetable.preference + phenofile_with_metab$score_group + phenofile_with_metab$sex)
summary(fit_hip_veg_pref)

fit_hip_fruit_pref <- lm(phenofile_with_metab$hippurate_log ~ phenofile_with_metab$fruit.preference)
summary(fit_hip_fruit_pref)
fit_hip_fruit_pref <- lm(phenofile_with_metab$hippurate_log ~ phenofile_with_metab$fruit.preference + phenofile_with_metab$sex)
summary(fit_hip_fruit_pref)
fit_hip_fruit_pref <- lm(phenofile_with_metab$hippurate_log ~ phenofile_with_metab$fruit.preference + phenofile_with_metab$score_group + phenofile_with_metab$sex)
summary(fit_hip_fruit_pref)

# pfos and fish
pfos_col <- which(names(phenofile_with_metab) == pfos_compid)
names(phenofile_with_metab)[pfos_col] <- "pfos"
# standardise
phenofile_with_metab$pfos_z <- scale(phenofile_with_metab$pfos)

fit_pfos_fish_pref <- lm(phenofile_with_metab$pfos_z ~ phenofile_with_metab$fish.preference)
summary(fit_pfos_fish_pref)
fit_pfos_fish_pref <- lm(phenofile_with_metab$pfos_z ~ phenofile_with_metab$fish.preference + phenofile_with_metab$sex)
summary(fit_pfos_fish_pref)
fit_pfos_fish_pref <- lm(phenofile_with_metab$pfos_z ~ phenofile_with_metab$fish.preference + phenofile_with_metab$score_group + phenofile_with_metab$sex)
summary(fit_pfos_fish_pref)
summary(fit_pfos_fish_pref)$coefficients[,4]

#################
###### END ######
#################

print(sessionInfo())

q()
