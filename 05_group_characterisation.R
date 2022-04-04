# script for analysing the association between adiposity traits and potential confounding factors and Yengo BMI GRS groups in RbG study
# also compared to an equivalent BMI phenotype recall study.

####################
###### SET UP ######
####################

# read in parameter file (specified on command line)
source("parameter_files/parameters_for_r.R")

# move to working directory 
setwd(working_dir)

## load libraries
library(heatmap3)
library(RColorBrewer)
library(tidyverse)
library(stringr)
library(data.table)
library(e1071)
library(ggsci)
library(ggplot2)

####################################################################################
####################################################################################
## read in data
####################################################################################
####################################################################################

## sample metadata with group added
orig_sample <- read.table(paste0(data_intermediate_dir,"Metabolon_sampledata_with_group.txt"), h=T, stringsAsFactors = F)

## pheno data ##
## extract from stata
phenofile <- read.table(file = paste0(data_intermediate_dir,"bmi_rbg_meta_processed_with_imputed_weight.csv"), h=T, sep=",", stringsAsFactors = F)                        

# read in GRS data used to do sample selection (RbG)
grs_data <- read.table(file="../selection01/Wade_20200317/yengo_score.profile",h=T)

####################################################################################
####################################################################################
## pheno process and checks
####################################################################################
####################################################################################

# check group allocation
table(orig_sample$score_group)

# add aln_qlet ID to phenofile (for use with merging)
phenofile$aln_qlet <- paste(phenofile$aln,phenofile$qlet,sep="_")
# add alnqlet to phenofile (for use with merging with score file)
phenofile$alnqlet <- paste(phenofile$aln,phenofile$qlet,sep="")

# convert negative data (various missing codes) to NAs in main phenofile (apply to variable columns only, not ID columns)
phenofile[,c(3:(dim(phenofile)[2]-2))][phenofile[,c(3:(dim(phenofile)[2]-2))] < 0] <- NA

# merge sample info with pheno, only keep those in metabolite sample file
combined_pheno_data <- merge(phenofile, orig_sample[,c("CLIENT_IDENTIFIER","score_group")], 
                             by.x="aln_qlet", by.y="CLIENT_IDENTIFIER", all.y="TRUE")
combined_pheno_data$score_group <- factor(combined_pheno_data$score_group)


####################################################################################
####################################################################################
## make data files for running models on
####################################################################################
####################################################################################

# check duplicates
duplicated_samples <- which(duplicated(combined_pheno_data$aln_qlet))
duplicated_samples

# drop duplicates that occur because samples from these individuals were analysed in duplicate
combined_pheno_data_nodups <-combined_pheno_data[-duplicated_samples,]
dim(combined_pheno_data_nodups)

# drop individuals removed during Metabolon metabolites QC
individuals_to_remove <- readRDS(file = paste0(data_intermediate_dir, "individuals_to_remove.rds"))
combined_pheno_data_nodups <- combined_pheno_data_nodups[!(combined_pheno_data_nodups$aln_qlet %in% individuals_to_remove),]
table(combined_pheno_data_nodups$score_group)

# clean blood pressure data
# derive measure to use as mean of the two measures taken at the given clinic (time point) for teenage clinics
combined_pheno_data_nodups$SBP_TF1 <- rowMeans(combined_pheno_data_nodups[,c("SBP1_TF1","SBP2_TF1")], na.rm = T)
combined_pheno_data_nodups$DBP_TF1 <- rowMeans(combined_pheno_data_nodups[,c("DBP1_TF1","DBP2_TF1")], na.rm = T)

combined_pheno_data_nodups$SBP_TF2 <- rowMeans(combined_pheno_data_nodups[,c("SBP1_TF2","SBP2_TF2")], na.rm = T)
combined_pheno_data_nodups$DBP_TF2 <- rowMeans(combined_pheno_data_nodups[,c("DBP1_TF2","DBP2_TF2")], na.rm = T)

combined_pheno_data_nodups$SBP_TF3 <- rowMeans(combined_pheno_data_nodups[,c("SBP1_TF3","SBP2_TF3")], na.rm = T)
combined_pheno_data_nodups$DBP_TF3 <- rowMeans(combined_pheno_data_nodups[,c("DBP1_TF3","DBP2_TF3")], na.rm = T)

# remove the measures 1 & 2
toMatch <- c("SBP1","SBP2","DBP1","DBP2")
combined_pheno_data_nodups <- combined_pheno_data_nodups[!grepl(paste(toMatch, collapse = "|"),names(combined_pheno_data_nodups))]
# replace NaN's generated above with NA
combined_pheno_data_nodups[combined_pheno_data_nodups=="NaN"] <- NA

clean_dat <- combined_pheno_data_nodups

write.table(clean_dat[!str_detect(names(clean_dat), "age_")],file = paste0(data_intermediate_dir, "phenotype_dat_clean.txt"),
            quote = F,row.names = F,sep="\t")
write.table(clean_dat,file = paste0(data_intermediate_dir, "phenotype_dat_clean_with_age.txt"),
            quote = F,row.names = F,sep="\t")


####################################################################################
####################################################################################
## extract stats to help with power calculations
####################################################################################
####################################################################################

### check GRS association with BMI
print("Relationship between GRS and BMI in entire eligible cohort with data")
grs_data$bmi <- phenofile$bmi_F24[match(grs_data$IID,phenofile$alnqlet)]
# remove missing BMI values
grs_data_complete <- grs_data[grs_data$bmi > 0,]
dim(grs_data_complete)
fit <- lm(bmi ~ SCORESUM, data=grs_data_complete)
summary(fit)
nobs(fit)

### check exposure-outcome associations using blood sample measures
fit <- lm(phenofile$Glucose_F24 ~ phenofile$bmi_F24)
summary(fit)
nobs(fit)

fit <- lm(phenofile$Chol_F24 ~ phenofile$bmi_F24)
summary(fit)
nobs(fit)

fit <- lm(phenofile$LDL_F24 ~ phenofile$bmi_F24)
summary(fit)
nobs(fit)

fit <- lm(phenofile$VLDL_F24 ~ phenofile$bmi_F24)
summary(fit)
nobs(fit)

fit <- lm(phenofile$Trig_F24 ~ phenofile$bmi_F24)
summary(fit)
nobs(fit)

fit <- lm(phenofile$Insulin_F24 ~ phenofile$bmi_F24)
summary(fit)
nobs(fit)

####################################################################################
####################################################################################
## sumarise basic info for pop charact table 
## Paper Ref: Table 1
####################################################################################
####################################################################################

# Note: 'phenofile' contains all F24 participants; 'clean_dat' contains only RbG participants.

### sex distribution (1=male, 2=female)
print("sex distribution (1=male, 2=female)")

# cohort at F24
print("Whole cohort")
table(phenofile$sex, useNA="always")

# in study
print("RbG study")
table(clean_dat$sex)

# high group
print("RbG study - high group")
table(clean_dat[clean_dat$score_group == "high", c("sex")])
# low group
print("RbG study - low group")
table(clean_dat[clean_dat$score_group == "low", c("sex")])

### fat mass
## transform to kg from g
phenofile$fat_mass_total_F24 <- phenofile$fat_mass_total_F24/1000
clean_dat$fat_mass_total_F24 <- clean_dat$fat_mass_total_F24/1000
### lean mass
## transform to kg from g
phenofile$lean_mass_total_F24 <- phenofile$lean_mass_total_F24/1000
clean_dat$lean_mass_total_F24 <- clean_dat$lean_mass_total_F24/1000
### age
## transform to years
phenofile$age_yrs_F24 <- phenofile$age_mth_F24/12
clean_dat$age_yrs_F24 <- clean_dat$age_mth_F24/12
### derive WHR
phenofile$whr_F24 <- phenofile$waist_circ_F24/phenofile$hip_circ_F24
clean_dat$whr_F24 <- clean_dat$waist_circ_F24/clean_dat$hip_circ_F24

# variables to summarise
char_var_list <- c("age_yrs_F24","bmi_F24","weight_F24","fat_mass_total_F24","lean_mass_total_F24","whr_F24")

# make results output
char_var_summ <- as.data.frame( matrix(data = NA, nrow = length(char_var_list), ncol = 25) )

# run analysis loop
for (i in 1:length(char_var_list)) {
  print(paste0("Processing summary for: ", char_var_list[i]))
  # var name
  char_var_summ[i,1] <- char_var_list[i]
  # whole cohort (phenofile)
  # var col
  col_num <- which(names(phenofile) == char_var_list[i])
  # whole cohort - all
  # calculate N 
  char_var_summ[i,2] <- length(which(!is.na(phenofile[,col_num])))
  # calculate means
  char_var_summ[i,3] <- mean(phenofile[,col_num], na.rm=T)
  # calculate SD
  char_var_summ[i,4] <- sd(phenofile[,col_num], na.rm=T)
  # whole cohort - males
  # calculate N 
  char_var_summ[i,5] <- length(which(!is.na(phenofile[phenofile$sex == 1,col_num])))
  # calculate means
  char_var_summ[i,6] <- mean(phenofile[phenofile$sex == 1,col_num], na.rm=T)
  # calculate SD
  char_var_summ[i,7] <- sd(phenofile[phenofile$sex == 1,col_num], na.rm=T)
  # whole cohort - females
  # calculate N 
  char_var_summ[i,8] <- length(which(!is.na(phenofile[phenofile$sex == 2,col_num])))
  # calculate means
  char_var_summ[i,9] <- mean(phenofile[phenofile$sex == 2,col_num], na.rm=T)
  # calculate SD
  char_var_summ[i,10] <- sd(phenofile[phenofile$sex == 2,col_num], na.rm=T)
  
  # RbG cohort (clean_dat)
  # var col
  col_num <- which(names(clean_dat) == char_var_list[i])
  # RbG cohort - all
  # calculate N 
  char_var_summ[i,11] <- length(which(!is.na(clean_dat[,col_num])))
  # calculate means
  char_var_summ[i,12] <- mean(clean_dat[,col_num], na.rm=T)
  # calculate SD
  char_var_summ[i,13] <- sd(clean_dat[,col_num], na.rm=T)
  # RbG cohort - males - high
  # calculate N 
  char_var_summ[i,14] <- length(which(!is.na(clean_dat[clean_dat$sex == 1 & clean_dat$score_group == "high",col_num])))
  # calculate means
  char_var_summ[i,15] <- mean(clean_dat[clean_dat$sex == 1 & clean_dat$score_group == "high",col_num], na.rm=T)
  # calculate SD
  char_var_summ[i,16] <- sd(clean_dat[clean_dat$sex == 1 & clean_dat$score_group == "high",col_num], na.rm=T)
  # RbG cohort - males - low
  # calculate N 
  char_var_summ[i,17] <- length(which(!is.na(clean_dat[clean_dat$sex == 1 & clean_dat$score_group == "low",col_num])))
  # calculate means
  char_var_summ[i,18] <- mean(clean_dat[clean_dat$sex == 1 & clean_dat$score_group == "low",col_num], na.rm=T)
  # calculate SD
  char_var_summ[i,19] <- sd(clean_dat[clean_dat$sex == 1 & clean_dat$score_group == "low",col_num], na.rm=T)
  # RbG cohort - females - high
  # calculate N 
  char_var_summ[i,20] <- length(which(!is.na(clean_dat[clean_dat$sex == 2 & clean_dat$score_group == "high",col_num])))
  # calculate means
  char_var_summ[i,21] <- mean(clean_dat[clean_dat$sex == 2 & clean_dat$score_group == "high",col_num], na.rm=T)
  # calculate SD
  char_var_summ[i,22] <- sd(clean_dat[clean_dat$sex == 2 & clean_dat$score_group == "high",col_num], na.rm=T)
  # RbG cohort - females - low
  # calculate N 
  char_var_summ[i,23] <- length(which(!is.na(clean_dat[clean_dat$sex == 2 & clean_dat$score_group == "low",col_num])))
  # calculate means
  char_var_summ[i,24] <- mean(clean_dat[clean_dat$sex == 2 & clean_dat$score_group == "low",col_num], na.rm=T)
  # calculate SD
  char_var_summ[i,25] <- sd(clean_dat[clean_dat$sex == 2 & clean_dat$score_group == "low",col_num], na.rm=T)
}

names(char_var_summ) <- c("variable","cohort_all_n","cohort_all_mean","cohort_all_sd","cohort_males_n","cohort_males_mean","cohort_males_sd","cohort_females_n","cohort_females_mean","cohort_females_sd",
                           "rbg_all_n","rbg_all_mean","rbg_all_sd",
                           "rbg_males_high_n","rbg_males_high_mean","rbg_males_high_sd","rbg_males_low_n","rbg_males_low_mean","rbg_males_low_sd",
                           "rbg_females_high_n","rbg_females_high_mean","rbg_females_high_sd","rbg_females_low_n","rbg_females_low_mean","rbg_females_low_sd")


# evaluate compliance
table(phenofile$fasting_status_f24, useNA = "always")
3103/(3103+360)

####################################################################################
####################################################################################
## summarise extra info for pop characteristics table 
## Paper Ref: Table S1
####################################################################################
####################################################################################

# Note: 'phenofile' contains all F24 participants; 'clean_dat' contains only RbG participants.

# variables to summarise
extra_var_list <- c("Insulin_F24","Glucose_F24","Chol_F24","Trig_F24","HDL_F24","LDL_F24","VLDL_F24","CRP_F24","seated_SBP_F24","seated_DBP_F24")

# make results output
extra_var_summ <- as.data.frame( matrix(data = NA, nrow = length(extra_var_list), ncol = 25) )

# run analysis loop
for (i in 1:length(extra_var_list)) {
  print(paste0("Processing summary for: ", extra_var_list[i]))
  # var name
  extra_var_summ[i,1] <- extra_var_list[i]
  # whole cohort (phenofile)
  # var col
  col_num <- which(names(phenofile) == extra_var_list[i])
  # whole cohort - all
  # calculate N 
  extra_var_summ[i,2] <- length(which(!is.na(phenofile[,col_num])))
  # calculate means
  extra_var_summ[i,3] <- mean(phenofile[,col_num], na.rm=T)
  # calculate SD
  extra_var_summ[i,4] <- sd(phenofile[,col_num], na.rm=T)
  # whole cohort - males
  # calculate N 
  extra_var_summ[i,5] <- length(which(!is.na(phenofile[phenofile$sex == 1,col_num])))
  # calculate means
  extra_var_summ[i,6] <- mean(phenofile[phenofile$sex == 1,col_num], na.rm=T)
  # calculate SD
  extra_var_summ[i,7] <- sd(phenofile[phenofile$sex == 1,col_num], na.rm=T)
  # whole cohort - females
  # calculate N 
  extra_var_summ[i,8] <- length(which(!is.na(phenofile[phenofile$sex == 2,col_num])))
  # calculate means
  extra_var_summ[i,9] <- mean(phenofile[phenofile$sex == 2,col_num], na.rm=T)
  # calculate SD
  extra_var_summ[i,10] <- sd(phenofile[phenofile$sex == 2,col_num], na.rm=T)
  
  # RbG cohort (clean_dat)
  # var col
  col_num <- which(names(clean_dat) == extra_var_list[i])
  # RbG cohort - all
  # calculate N 
  extra_var_summ[i,11] <- length(which(!is.na(clean_dat[,col_num])))
  # calculate means
  extra_var_summ[i,12] <- mean(clean_dat[,col_num], na.rm=T)
  # calculate SD
  extra_var_summ[i,13] <- sd(clean_dat[,col_num], na.rm=T)
  # RbG cohort - males - high
  # calculate N 
  extra_var_summ[i,14] <- length(which(!is.na(clean_dat[clean_dat$sex == 1 & clean_dat$score_group == "high",col_num])))
  # calculate means
  extra_var_summ[i,15] <- mean(clean_dat[clean_dat$sex == 1 & clean_dat$score_group == "high",col_num], na.rm=T)
  # calculate SD
  extra_var_summ[i,16] <- sd(clean_dat[clean_dat$sex == 1 & clean_dat$score_group == "high",col_num], na.rm=T)
  # RbG cohort - males - low
  # calculate N 
  extra_var_summ[i,17] <- length(which(!is.na(clean_dat[clean_dat$sex == 1 & clean_dat$score_group == "low",col_num])))
  # calculate means
  extra_var_summ[i,18] <- mean(clean_dat[clean_dat$sex == 1 & clean_dat$score_group == "low",col_num], na.rm=T)
  # calculate SD
  extra_var_summ[i,19] <- sd(clean_dat[clean_dat$sex == 1 & clean_dat$score_group == "low",col_num], na.rm=T)
  # RbG cohort - females - high
  # calculate N 
  extra_var_summ[i,20] <- length(which(!is.na(clean_dat[clean_dat$sex == 2 & clean_dat$score_group == "high",col_num])))
  # calculate means
  extra_var_summ[i,21] <- mean(clean_dat[clean_dat$sex == 2 & clean_dat$score_group == "high",col_num], na.rm=T)
  # calculate SD
  extra_var_summ[i,22] <- sd(clean_dat[clean_dat$sex == 2 & clean_dat$score_group == "high",col_num], na.rm=T)
  # RbG cohort - females - low
  # calculate N 
  extra_var_summ[i,23] <- length(which(!is.na(clean_dat[clean_dat$sex == 2 & clean_dat$score_group == "low",col_num])))
  # calculate means
  extra_var_summ[i,24] <- mean(clean_dat[clean_dat$sex == 2 & clean_dat$score_group == "low",col_num], na.rm=T)
  # calculate SD
  extra_var_summ[i,25] <- sd(clean_dat[clean_dat$sex == 2 & clean_dat$score_group == "low",col_num], na.rm=T)
}

names(extra_var_summ) <- c("variable","cohort_all_n","cohort_all_mean","cohort_all_sd","cohort_males_n","cohort_males_mean","cohort_males_sd","cohort_females_n","cohort_females_mean","cohort_females_sd",
                           "rbg_all_n","rbg_all_mean","rbg_all_sd",
                           "rbg_males_high_n","rbg_males_high_mean","rbg_males_high_sd","rbg_males_low_n","rbg_males_low_mean","rbg_males_low_sd",
                           "rbg_females_high_n","rbg_females_high_mean","rbg_females_high_sd","rbg_females_low_n","rbg_females_low_mean","rbg_females_low_sd")


####################################################################################
####################################################################################
## assess potential covariates
####################################################################################
####################################################################################

conti_covar <- c("MVPA_minutes_F24","age_wks_F24","time_since_last_food_F24")
nonconti_covar <- c("sex","mum_ever_smoke","mum_alc_freq","mum_highest_edu",
                    "mat_social_class","pat_social_class","parity","fasting_status_f24")
covar_dat <- clean_dat[,c(conti_covar, nonconti_covar, "score_group")]
covar_dat[conti_covar] <- sapply(covar_dat[conti_covar], as.numeric)
covar_dat[nonconti_covar] <- lapply(covar_dat[nonconti_covar], factor)

binary_covar <- c("sex","mum_ever_smoke","fasting_status_f24")
categorical_covar <- c("mum_highest_edu","mat_social_class","pat_social_class","parity","mum_alc_freq")

## perform tests
# continous variable
conti_res <- data.frame()
for (i in 1:length(conti_covar)) {
  var_name <- conti_covar[i]
  df_tmp <- na.omit(covar_dat[,c(var_name,"score_group")])
  if (levels(df_tmp$score_group)[1] == "low") {
    df_tmp <- within(df_tmp,score_group <- relevel(score_group, ref = "high"))
  }
  shapwilk <- shapiro.test(df_tmp[,var_name])
  res <- t.test(as.formula(paste0(var_name,"~","score_group")),data = df_tmp)
  conti_res[i,1:6] <- c(var_name, 
                            res$estimate[1] - res$estimate[2], 
                            res$p.value,
                            res$conf.int[1], 
                            res$conf.int[2], 
                            res$parameter)
  conti_res[i,7:8] <- c(shapwilk$statistic,shapwilk$p.value) 
}
colnames(conti_res) <- c("phenotype","mean_difference","pval","ci_low","ci_high","df","shapirowilk_W","shapirowilk_p")

# binary variable
binary_res <- data.frame()
for (i in 1:length(binary_covar)) {
  var_name <- binary_covar[i]
  print(paste0(var_name))
  print(table(covar_dat[,c(var_name)]))
  df_tmp <- na.omit(covar_dat[,c(var_name,"score_group")])
  if (levels(df_tmp$score_group)[1] == "low") {
    df_tmp <- within(df_tmp,score_group <- relevel(score_group, ref = "high"))
  }
  res <- fisher.test(table(df_tmp[[var_name]],df_tmp$score_group))
  binary_res[i,1:6] <- c(var_name, 
                         res$p.value, 
                         dim(df_tmp)[1], 
                         res$estimate, 
                         res$conf.int[1], 
                         res$conf.int[2])
}
names(binary_res) <- c("phenotype", "pval", "N", "OR", "OR_CI_low", "OR_CI_high")

# categorical variables
categorical_res <- data.frame()
for (i in 1:length(categorical_covar)) {
  var_name <- categorical_covar[i]
  df_tmp <- na.omit(covar_dat[,c(var_name,"score_group")])
  if (levels(df_tmp$score_group)[1] == "low") {
    df_tmp <- within(df_tmp,score_group <- relevel(score_group, ref = "high"))
  }
  res <- fisher.test(table(df_tmp[[var_name]],df_tmp$score_group),simulate.p.value = TRUE)
  categorical_res[i,1:3] <- c(var_name, 
                              res$p.value, 
                              dim(df_tmp)[1])
}
names(categorical_res) <- c("phenotype", "pval", "N")

# Paper Ref: Table S2 [manually compiled from the following]
fwrite(conti_res, paste0(data_output_dir, "ForPaper/TableS3_05.1_confounder_analysis_continious_variables.txt"), sep = "\t")
fwrite(binary_res, paste0(data_output_dir, "ForPaper/TableS3_05.1_confounder_analysis_binary_variables.txt"), sep = "\t")
fwrite(categorical_res, paste0(data_output_dir, "ForPaper/TableS3_05.1_confounder_analysis_categorical_variables.txt"), sep = "\t")

## create tables for variables with p<0.10
## cross tab mat and pat social class with score group
mat_social_class_temp <- table(covar_dat$mat_social_class,covar_dat$score_group)
pat_social_class_temp <- table(covar_dat$pat_social_class,covar_dat$score_group)

# put table counts data in dataframe
social_class_tab <- data.frame(c("professional", "managerial/technical", "skilled (non-manual", "skilled (manual)", "partly-skilled", "unskilled", "armed forces"), 
                               c(mat_social_class_temp[1:6],"0"), c(mat_social_class_temp[7:12],"0"), pat_social_class_temp[1:7], pat_social_class_temp[8:14])

# rename columns
names(social_class_tab)[1] <- "social_class"
names(social_class_tab)[2] <- "mat_high"
names(social_class_tab)[3] <- "mat_low"
names(social_class_tab)[4] <- "pat_high"
names(social_class_tab)[5] <- "pat_low"

# calculate %
social_class_tab$mat_high_pct <- round(as.numeric(social_class_tab$mat_high)/sum(as.numeric(social_class_tab$mat_high)), digits=3)
social_class_tab$mat_low_pct <- round(as.numeric(social_class_tab$mat_low)/sum(as.numeric(social_class_tab$mat_low)), digits=3)
social_class_tab$pat_high_pct <- round(as.numeric(social_class_tab$pat_high)/sum(as.numeric(social_class_tab$pat_high)), digits=3)
social_class_tab$pat_low_pct <- round(as.numeric(social_class_tab$pat_low)/sum(as.numeric(social_class_tab$pat_low)), digits=3)

# save out
fwrite(social_class_tab, paste0(data_output_dir, "05.1_confounder_analysis_social_class_table.txt"), sep = "\t")

# graph social class data
graph_data_class <- as.vector(social_class_tab[,c("social_class")])
# maternal
graph_data_mat_low <- as.vector(social_class_tab[,c("mat_low")])
graph_data_mat_high <- as.vector(social_class_tab[,c("mat_high")])
graph_data <- data.frame(class=c(graph_data_class, graph_data_class),count=c(graph_data_mat_low, graph_data_mat_high),
                         score_group = c("low.BMI.GRS","low.BMI.GRS","low.BMI.GRS","low.BMI.GRS","low.BMI.GRS","low.BMI.GRS","low.BMI.GRS",
                                         "high.BMI.GRS","high.BMI.GRS","high.BMI.GRS","high.BMI.GRS","high.BMI.GRS","high.BMI.GRS","high.BMI.GRS"))
graph_data$score_group <- as.factor(graph_data$score_group)
graph_data$class <- as.factor(graph_data$class)
graph_data$count <- as.numeric(graph_data$count)


ggplot(graph_data, 
       aes(x=factor(class, levels = c("armed forces", "unskilled", "partly-skilled", "skilled (manual)", "skilled (non-manual", "managerial/technical","professional" )),
           y= count, fill = factor(score_group))) + 
  geom_bar(stat = "identity", position = position_dodge()) + 
  theme(axis.text.x = element_text(angle = 45, 
                                   hjust = 1)) +
  labs(x = "Maternal social class", y= "Frequency", fill = "BMI GRS group")
ggsave(paste0(data_output_dir, "ForPaper/FigS2a_05.1_mat_social_class.pdf"), 
       plot = last_plot(),
       height = 10, width = 10, units = "cm",
       dpi = 320)

# paternal
graph_data_pat_low <- as.vector(social_class_tab[,c("pat_low")])
graph_data_pat_high <- as.vector(social_class_tab[,c("pat_high")])
graph_data <- data.frame(class=c(graph_data_class, graph_data_class),count=c(graph_data_pat_low, graph_data_pat_high),
                         score_group = c("low.BMI.GRS","low.BMI.GRS","low.BMI.GRS","low.BMI.GRS","low.BMI.GRS","low.BMI.GRS","low.BMI.GRS",
                                         "high.BMI.GRS","high.BMI.GRS","high.BMI.GRS","high.BMI.GRS","high.BMI.GRS","high.BMI.GRS","high.BMI.GRS"))
graph_data$score_group <- as.factor(graph_data$score_group)
graph_data$class <- as.factor(graph_data$class)
graph_data$count <- as.numeric(graph_data$count)


ggplot(graph_data, 
       aes(x=factor(class, levels = c("armed forces", "unskilled", "partly-skilled", "skilled (manual)", "skilled (non-manual", "managerial/technical","professional" )),
           y= count, fill = factor(score_group))) + 
  geom_bar(stat = "identity", position = position_dodge()) + 
  theme(axis.text.x = element_text(angle = 45, 
                                   hjust = 1)) +
  labs(x = "Paternal social class", y= "Frequency", fill = "BMI GRS group")
ggsave(paste0(data_output_dir, "ForPaper/FigS2b_05.1_pat_social_class.pdf"), 
       plot = last_plot(),
       height = 10, width = 10, units = "cm",
       dpi = 320)


####################################################################################
####################################################################################
## create age and batch vectors
####################################################################################
####################################################################################

age_dat <- combined_pheno_data_nodups[str_detect(names(combined_pheno_data_nodups), "age_wks")]
yrs_to_wks <- 52.17857
mth_to_wks <- 4.34524

age_means <- colMeans(age_dat, na.rm = T)
age_means <- unlist(age_means)
age_means <- age_means[c(9:18,8:1,19:20)]

names(age_means) <- gsub("age_wks_","",names(age_means))

batch_list <- c(rep("Child in Focus",10),rep("Focus@7-11",5),rep("Teen Focus",4),"Focus@24")
names(batch_list) <- gsub("age_wks_","",names(age_means))

####################################################################################
####################################################################################
## assessing positive control phenotypes (BMI, weight, blood pressure etc.)
####################################################################################
####################################################################################

pheno_list <- c("age_yrs_F24","weight","bmi","birthweight","waist_circ","fat_mass_total","lean_mass_total","whr_F24",
                "Glucose","Chol","Trig","HDL","LDL","VLDL","CRP","Insulin","SBP","DBP")
var_list <- names(clean_dat)[grepl(paste(c(pheno_list), collapse = "|"), names(clean_dat))]

dtst <- clean_dat[c(var_list, "score_group")]

# perform t-test on phenotype measures
res_mean_difference <- data.frame()
for (i in 1:length(var_list)) {
  var_name <- var_list[i]
  res_mean_difference[i,1] <- var_name
  df_tmp <- na.omit(dtst[c(var_name, "score_group")])
  if (levels(df_tmp$score_group)[1] == "low") {
    df_tmp <- within(df_tmp,score_group <- relevel(score_group, ref = "high"))
  }
  # check distribution
  shapwilk <- shapiro.test(df_tmp[,var_name])
  res <- t.test(as.formula(paste(var_name, "score_group", sep = "~")),
                data = df_tmp)
  res_mean_difference[i,2:6] <- c(res$estimate[1] - res$estimate[2],
                                  res$p.value,
                                  res$conf.int[1],
                                  res$conf.int[2],
                                  res$parameter)
  res_mean_difference[i,7:8] <- c(shapwilk$statistic,shapwilk$p.value) 
  res_mean_difference[i,9] <- nrow(df_tmp)
  # add in rank test
  res_mean_difference[i,10] <- wilcox.test(as.formula(paste(var_name, "score_group", sep = "~")),
                                           data = df_tmp)$p.value
}
colnames(res_mean_difference) <- c("phenotype","mean_difference","pval","ci_low","ci_high","df",
                                   "shapirowilk_W","shapirowilk_p","N","wilcox_pval")

# re-order data by age
res_mean_difference <- rbind(res_mean_difference[res_mean_difference$phenotype == "birthweight",],
                             res_mean_difference[grepl("mth",res_mean_difference$phenotype), ], 
                             res_mean_difference[grepl("F7",res_mean_difference$phenotype), ], 
                             res_mean_difference[grepl("F8",res_mean_difference$phenotype), ], 
                             res_mean_difference[grepl("F9",res_mean_difference$phenotype), ], 
                             res_mean_difference[grepl("F10",res_mean_difference$phenotype), ], 
                             res_mean_difference[grepl("F11",res_mean_difference$phenotype), ], 
                             res_mean_difference[grepl("TF1",res_mean_difference$phenotype), ], 
                             res_mean_difference[grepl("TF2",res_mean_difference$phenotype), ], 
                             res_mean_difference[grepl("TF3",res_mean_difference$phenotype), ], 
                             res_mean_difference[grepl("F17",res_mean_difference$phenotype), ], 
                             res_mean_difference[grepl("F24",res_mean_difference$phenotype), ])

fwrite(res_mean_difference, paste0(data_output_dir, "05.2_phenotype_analysis_mean_differences_results.txt"), sep = "\t")

## join test stats to summary tables from above and write out
char_var_summ_plus <- left_join(char_var_summ,res_mean_difference,by = c("variable" = "phenotype"))
extra_var_summ_plus <- left_join(extra_var_summ,res_mean_difference,by = c("variable" = "phenotype"))

# save out info for Table 1
fwrite(char_var_summ_plus, paste0(data_output_dir, "ForPaper/Table1_05.3_characterisation.txt"), sep = "\t")

# update variable names to be more reader friendly
extra_var_summ_plus$variable
extra_var_summ_plus <- extra_var_summ_plus %>% mutate(`variable` = replace(`variable`, `variable` == "Insulin_F24", "Fasting.insulin"))
extra_var_summ_plus <- extra_var_summ_plus %>% mutate(`variable` = replace(`variable`, `variable` == "Glucose_F24", "Fasting.glucose"))
extra_var_summ_plus <- extra_var_summ_plus %>% mutate(`variable` = replace(`variable`, `variable` == "Chol_F24", "Fasting.cholesterol"))
extra_var_summ_plus <- extra_var_summ_plus %>% mutate(`variable` = replace(`variable`, `variable` == "LDL_F24", "LDL-C"))
extra_var_summ_plus <- extra_var_summ_plus %>% mutate(`variable` = replace(`variable`, `variable` == "VLDL_F24", "VLDL-C"))
extra_var_summ_plus <- extra_var_summ_plus %>% mutate(`variable` = replace(`variable`, `variable` == "HDL_F24", "HDL-C"))
extra_var_summ_plus <- extra_var_summ_plus %>% mutate(`variable` = replace(`variable`, `variable` == "CRP_F24", "C-reactive.protein"))
extra_var_summ_plus <- extra_var_summ_plus %>% mutate(`variable` = replace(`variable`, `variable` == "Trig_F24", "Triglycerides"))
extra_var_summ_plus <- extra_var_summ_plus %>% mutate(`variable` = replace(`variable`, `variable` == "seated_SBP_F24", "Seated.systolic.blood.pressure"))
extra_var_summ_plus <- extra_var_summ_plus %>% mutate(`variable` = replace(`variable`, `variable` == "seated_DBP_F24", "Seated.diastolic.blood.pressure"))
extra_var_summ_plus$variable

# save out table S1
fwrite(extra_var_summ_plus, paste0(data_output_dir, "ForPaper/TableS1_05.3_extra_characterisation.txt"), sep = "\t")


####################################################################################
####################################################################################
## extract weight and bmi temporal data for suppl table
####################################################################################
####################################################################################

# Paper Ref: Table S2

# extract rows for bmi and weight
bmi_rows <- res_mean_difference[grep("bmi", res_mean_difference$phenotype),]
weight_rows <- res_mean_difference[grep("weight_", res_mean_difference$phenotype),]
bmi_weight_rows <- rbind(bmi_rows,weight_rows)
# separate clinic label
bmi_weight_out <- bmi_weight_rows %>% separate(phenotype,c("measure","clinic"),sep="_")
# add mean age for each clinic
age_means_dta <- as.data.frame(age_means)
age_means_dta$clinic_identifier <- rownames(age_means_dta)
bmi_weight_out$mean_age_in_weeks <- age_means_dta$age_means[match(bmi_weight_out$clinic,age_means_dta$clinic)]
# transform to years
bmi_weight_out$mean_age_in_years <- bmi_weight_out$mean_age_in_weeks/52

# re-order
bmi_weight_out <- bmi_weight_out[,c(1,2,11,12,10,3,5,6,4,7,8,9)]
# tidy up column headers
names(bmi_weight_out) <- c("Adiposity.measure", "ALSPAC.clinic.identifier","Mean.age.(weeks)","Mean.age.(years)","Sample.size",
                            "Between.group.difference.in.means.(high.cf.low)","Lower.95%.CI","Upper.95%.CI","Ttest.p.value","Degrees.of.freedom",
                            "Shapiro.Wilk.Wstat","Shapiro.Wilk.p.value")

# write out
names(bmi_weight_out)
fwrite(bmi_weight_out, paste0(data_output_dir, "ForPaper/TableS2_05.2_mean_differences_results.txt"), sep = "\t")

# plotting mean differences of variables across time
var_list <- c("weight","bmi")
unit_list <- c("kg", "kg/m2")

for (i in 1:length(var_list)) {
  var_name <- var_list[i]
  unit_name <- unit_list[i]
  df <- res_mean_difference %>% filter(str_detect(phenotype,var_name))
  df$time <- age_means[gsub(paste(var_name,"_",sep = ""),"",df$phenotype)] / yrs_to_wks
  df$pval_cat <- ifelse(df$pval < 0.001, "***", ifelse(df$pval < 0.01, "**", ifelse(df$pval < 0.1, "*", NA)))
  distance <- (range(df$ci_high)[2]-range(df$ci_low)[1]) * 0.02
  # get bmi in upper case
  var_label = NULL
  if(var_name == "bmi") {var_label = "BMI"}
  if(var_name == "weight") {var_label = "weight"}
  ggplot(data=df, aes(x=time, y=mean_difference, label=pval_cat)) +
    geom_point() +
    geom_text(hjust = 0.5, size=4, aes(y = ci_high + distance), show.legend = FALSE) +
    geom_errorbar(aes(ymin=ci_low,ymax=ci_high),width=0.5) +
    scale_x_continuous(breaks = seq(0,25,5), minor_breaks = seq(0 , 25, 1)) +
    scale_y_continuous(breaks = seq(0,15,2.5), minor_breaks = NULL) +
    geom_hline(yintercept = 0, linetype="dashed", color = "dark gray") +
    labs(color="Age groups", x="Age at clinic visit (years)", y= paste0("Mean differences in ", var_label, " (", unit_name,")")) +
    theme_bw() + scale_color_jco() + scale_fill_jco()
  ggsave(paste0(data_output_dir, "ForPaper/Fig3_and_S1_05.3_RbG_mean_difference_",var_name,".pdf"), 
         plot = last_plot(),
         height = 10, width = 20, units = "cm",
         dpi = 320)
}

####################################################################################
####################################################################################
## assessing confounding effects in recall-by-phenotype (BMI) result
####################################################################################
####################################################################################

# read in data
bmi_sample <- read.table(paste0(data_intermediate_dir,"BMI_sampledata_with_group.txt"), h=T, stringsAsFactors = F)

# check group allocation
table(bmi_sample$bmi_group)


combined_pheno_data <- merge(phenofile, bmi_sample[,c("aln_qlet","bmi_group")], 
                             by ="aln_qlet", all.y="TRUE")
combined_pheno_data$bmi_group <- factor(combined_pheno_data$bmi_group)

# check duplicates
duplicated_samples <- which(duplicated(combined_pheno_data$aln_qlet))
duplicated_samples

# drop duplicates
combined_pheno_data_nodups <- combined_pheno_data
dim(combined_pheno_data_nodups)
table(combined_pheno_data_nodups$bmi_group)

# convert negative data (missing) to NAs
combined_pheno_data_nodups[,c(4:(dim(combined_pheno_data_nodups)[2]-1))][combined_pheno_data_nodups[,c(4:(dim(combined_pheno_data_nodups)[2]-1))] < 0] <- NA

clean_dat <- combined_pheno_data_nodups

conti_covar <- c("MVPA_minutes_F24","age_wks_F24","time_since_last_food_F24")
nonconti_covar <- c("sex","mum_ever_smoke","mum_alc_freq","mum_highest_edu",
                    "mat_social_class","pat_social_class","parity","fasting_status_f24")
covar_dat <- clean_dat[,c(conti_covar, nonconti_covar, "bmi_group")]
covar_dat[conti_covar] <- sapply(covar_dat[conti_covar], as.numeric)
covar_dat[nonconti_covar] <- lapply(covar_dat[nonconti_covar], factor)

binary_covar <- c("sex","mum_ever_smoke","fasting_status_f24")
categorical_covar <- c("mum_highest_edu","mat_social_class","pat_social_class","parity","mum_alc_freq")

## perform tests
# continuous variable
conti_res <- data.frame()
for (i in 1:length(conti_covar)) {
  var_name <- conti_covar[i]
  df_tmp <- na.omit(covar_dat[,c(var_name,"bmi_group")])
  if (levels(df_tmp$bmi_group)[1] == "low") {
    df_tmp <- within(df_tmp,bmi_group <- relevel(bmi_group, ref = "high"))
  }
  shapwilk <- shapiro.test(df_tmp[,var_name])
  res <- t.test(as.formula(paste0(var_name,"~","bmi_group")),data = df_tmp)
  conti_res[i,1:6] <- c(var_name, 
                        res$estimate[1] - res$estimate[2], 
                        res$p.value,
                        res$conf.int[1], 
                        res$conf.int[2], 
                        res$parameter)
  conti_res[i,7:8] <- c(shapwilk$statistic,shapwilk$p.value) 
}
colnames(conti_res) <- c("phenotype","mean_difference","pval","ci_low","ci_high","df","shapirowilk_W","shapirowilk_p")

# binary variable
binary_res <- data.frame()
for (i in 1:length(binary_covar)) {
  var_name <- binary_covar[i]
  df_tmp <- na.omit(covar_dat[,c(var_name,"bmi_group")])
  if (levels(df_tmp$bmi_group)[1] == "low") {
    df_tmp <- within(df_tmp,bmi_group <- relevel(bmi_group, ref = "high"))
  }
  res <- fisher.test(table(df_tmp[[var_name]],df_tmp$bmi_group))
  binary_res[i,1:6] <- c(var_name, 
                         res$p.value, 
                         dim(df_tmp)[1], 
                         res$estimate, 
                         res$conf.int[1], 
                         res$conf.int[2])
}
names(binary_res) <- c("phenotype", "pval", "N", "OR", "OR_CI_low", "OR_CI_high")

# categorical variables
categorical_res <- data.frame()
for (i in 1:length(categorical_covar)) {
  var_name <- categorical_covar[i]
  df_tmp <- na.omit(covar_dat[,c(var_name,"bmi_group")])
  if (levels(df_tmp$bmi_group)[1] == "low") {
    df_tmp <- within(df_tmp,bmi_group <- relevel(bmi_group, ref = "high"))
  }
  res <- fisher.test(table(df_tmp[[var_name]],df_tmp$bmi_group),simulate.p.value = TRUE)
  categorical_res[i,1:3] <- c(var_name, 
                              res$p.value, 
                              dim(df_tmp)[1])
}
names(categorical_res) <- c("phenotype", "pval", "N")

fwrite(conti_res, paste0(data_output_dir, "05.4_RbBMI_confounder_analysis_continious_variables.txt"), sep = "\t")
fwrite(binary_res, paste0(data_output_dir, "05.4_RbBMI_confounder_analysis_binary_variables.txt"), sep = "\t")
fwrite(categorical_res, paste0(data_output_dir, "05.4_RbBMI_confounder_analysis_categorical_variables.txt"), sep = "\t")


## create tables for variables with p<0.10
## cross tab mat and pat social class with score group
mat_social_class_temp <- table(covar_dat$mat_social_class,covar_dat$bmi_group)
pat_social_class_temp <- table(covar_dat$pat_social_class,covar_dat$bmi_group)
mums_ed_temp <- table(covar_dat$mum_highest_edu,covar_dat$bmi_group)

# put table counts data in dataframes
social_class_tab <- data.frame(c("professional", "managerial/technical", "skilled (non-manual", "skilled (manual)", "partly-skilled", "unskilled", "armed forces"), c(mat_social_class_temp[1:6],"0"), c(mat_social_class_temp[7:12],"0"),
                               pat_social_class_temp[1:7], pat_social_class_temp[8:14])

mums_ed_tab <- data.frame(c("CSE", "Vocational", "O level", "A level", "Degree"), mums_ed_temp[1:5], mums_ed_temp[6:10])
# rename columns
names(social_class_tab)[1] <- "social_class"
names(social_class_tab)[2] <- "mat_high"
names(social_class_tab)[3] <- "mat_low"
names(social_class_tab)[4] <- "pat_high"
names(social_class_tab)[5] <- "pat_low"

names(mums_ed_tab)[1] <- "mum_highest_edu"
names(mums_ed_tab)[2] <- "high"
names(mums_ed_tab)[3] <- "low"

# calculate %
social_class_tab$mat_high_pct <- round(as.numeric(social_class_tab$mat_high)/sum(as.numeric(social_class_tab$mat_high)), digits=3)
social_class_tab$mat_low_pct <- round(as.numeric(social_class_tab$mat_low)/sum(as.numeric(social_class_tab$mat_low)), digits=3)
social_class_tab$pat_high_pct <- round(as.numeric(social_class_tab$pat_high)/sum(as.numeric(social_class_tab$pat_high)), digits=3)
social_class_tab$pat_low_pct <- round(as.numeric(social_class_tab$pat_low)/sum(as.numeric(social_class_tab$pat_low)), digits=3)

mums_ed_tab$high_pct <- round(as.numeric(mums_ed_tab$high)/sum(as.numeric(mums_ed_tab$high)), digits=3)
mums_ed_tab$low_pct <- round(as.numeric(mums_ed_tab$low)/sum(as.numeric(mums_ed_tab$low)), digits=3)

# save out
fwrite(social_class_tab, paste0(data_output_dir, "05.4_RbBMI_confounder_analysis_social_class_table.txt"), sep = "\t")
fwrite(mums_ed_tab, paste0(data_output_dir, "05.4_RbBMI_confounder_analysis_mums_ed_table.txt"), sep = "\t")

####################################################################################

#################
###### END ######
#################

print(sessionInfo())

q()
