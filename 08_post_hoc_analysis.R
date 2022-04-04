# script for post hoc analysis

####################
###### SET UP ######
####################

source("00_functions.R")

# read in parameter file (specified on command line)
source("parameter_files/parameters_for_r.R")

# move to working directory 
setwd(working_dir)

## load libraries
library(car)
library(heatmap3)
library(RColorBrewer)
library(tidyverse)
library(patchwork)
library(data.table)
library(iPVs)
library(ggpubr)
library(broom)
library(ggsci)
library(lmtest)
library(sandwich)
library(reshape2)
library(DescTools)
library(ggplot2)
library(dplyr)

####################################################################################
####################################################################################
## read in data
####################################################################################
####################################################################################

## Metabolon data ##
## original scale version (after QC)
orig_data <- read.table(paste0(metab_input_dir,"qc_data/B3194_BMI_RbG_2020_02_03_QCd_data.txt"), h=T)

## read updated feature metadata (new IDs and corrected labels included)
orig_feature <- read.table(paste0(data_intermediate_dir,"06.2_updated_orig_feature.txt"), 
                           h=T, stringsAsFactors = F,sep="\t")
orig_feature$feature_id <- rownames(orig_feature)

## sample metadata with group added
orig_sample <- read.table(paste0(data_intermediate_dir,"Metabolon_sampledata_with_group.txt"), 
                          h=T, stringsAsFactors = F)

## pheno data ##
## extract from 05 script
phenofile <- read.table(file = paste0(data_intermediate_dir,"phenotype_dat_clean_with_age.txt"), 
                        h=T, stringsAsFactors = F)

## linear regression results
merged_res <- fread(paste0(data_output_dir,"06.4_linear_analysis_results_all_models.txt"), sep = "\t")

## logistic regression results
logistic_res <- fread(paste0(data_output_dir,"07.1_logistic_regression_results_simple_model.txt"), sep = "\t")

## read in list of ALSPAC IDs to remove because they're part of a quad set
#quad_exc_list <- read.table(paste0(data_intermediate_dir, "quad_exc_list.txt"), h=F, stringsAsFactors = F)

## read in list of ALSPAC IDs to remove 
individuals_to_remove <- readRDS(file = paste0(data_intermediate_dir, "individuals_to_remove.rds"))

####################################################################################
####################################################################################
## pheno process and checks
####################################################################################
####################################################################################

# check group allocation
table(orig_sample$score_group)

# add aln_qlet ID to phenofile (for use with merging)
phenofile$aln_qlet <- paste(phenofile$aln,phenofile$qlet,sep="_")
head(phenofile)

# remove score_group column
phenofile <- subset(phenofile, select = -c(score_group))

####################################################################################
####################################################################################
## make data files for running models on
####################################################################################
####################################################################################

## generate feature list
feature_list <- names(orig_data)[1:(ncol(orig_data))]
head(feature_list)

# bring in sample id
orig_data$SAMPLE_NAME <- rownames(orig_data)

# drop duplicates
sample_to_remove <- readRDS(file = paste0(data_intermediate_dir, "samples_to_remove.rds"))
metab_nodups <- orig_sample[-which(orig_sample$SAMPLE_NAME %in% sample_to_remove),]
dim(metab_nodups)
table(metab_nodups$score_group)

# merge, keeping only those that are in both files
metab_full <- merge(metab_nodups,orig_data,by="SAMPLE_NAME")
dim(metab_full)

# change group to factor
metab_full$score_group <- as.factor(metab_full$score_group)
rownames(metab_full) <- metab_full$CLIENT_IDENTIFIER

# drop individuals to remove
metab_full <- metab_full[-which(metab_full$CLIENT_IDENTIFIER %in% individuals_to_remove),]
table(metab_full$score_group)

# merge sample info with pheno, only keep those in metabolite sample file
combined_pheno_data <- merge(phenofile, metab_nodups, by.x="aln_qlet", 
                             by.y="CLIENT_IDENTIFIER", all=FALSE)

## extract bmi data
bmi_data <- combined_pheno_data[,c("aln_qlet","bmi_F24","score_group")]

## extract lists of associated metabolites from main analyses
metab_list <- unlist(merged_res[which(merged_res$Associated.flag==1),]$feature_id)
# write out to use in cross validation
write.table(metab_list, paste0(data_output_dir, "08.2_linear_reg_assoc_list.txt"),row.names=F,quote=F,col.names=F)
metab_list_logistic <- unlist(logistic_res[which(logistic_res$logistic_flag==1),]$feature_id)

## restrict metabolite dataset to those associated metabolites from linear model
dtst <- merge(metab_full[,c("CLIENT_IDENTIFIER",metab_list)], 
              combined_pheno_data[,c("aln_qlet","fat_mass_total_F24","lean_mass_total_F24",
                                     "bmi_F24","score_group","sex","age_wks_F24")],
              by.x = "CLIENT_IDENTIFIER", by.y = "aln_qlet", all = FALSE)
# scale and mean centre
dtst <- dtst %>% mutate_at(metab_list, ~(scale(.) %>% as.vector))
dtst$score_group <- factor(dtst$score_group, levels = c("low","high"))
dtst$sex <- as.factor(dtst$sex)

rownames(dtst) <- dtst$CLIENT_IDENTIFIER

## restrict metabolite dataset to those associated metabolites from linear model and logistic model
dtst_with_logistic <- merge(metab_full[,c("CLIENT_IDENTIFIER",metab_list,metab_list_logistic)], 
                            combined_pheno_data[,c("aln_qlet","fat_mass_total_F24","lean_mass_total_F24",
                                                   "bmi_F24","score_group","sex","age_wks_F24")],
                            by.x = "CLIENT_IDENTIFIER", by.y = "aln_qlet", all = FALSE)
# scale and mean centre
dtst_with_logistic <- dtst_with_logistic %>% mutate_at(c(metab_list,metab_list_logistic), ~(scale(.) %>% as.vector))
dtst_with_logistic$score_group <- factor(dtst_with_logistic$score_group, levels = c("low","high"))
rownames(dtst_with_logistic) <- dtst_with_logistic$CLIENT_IDENTIFIER

####################################################################################
####################################################################################
## check data distributions by score group
####################################################################################
####################################################################################

# Box and whisker plots for all metabolites from linear regression results #
# Based on unimputed, Z-scored metabolite data

# Paper Ref: Fig S3

pdf_out = paste0(data_output_dir,"ForPaper/FigS3_08.0_boxplots_of_metabolite_levels_by_group.pdf")
pdf(pdf_out)
par(mfrow=c(4,3))
par(mar = c(2,5,2,2))

## plot
for (i in 2:(ncol(dtst)-6) ){
  metabolite_id <- names(dtst)[i]
  metabolite_name <- ifelse(metabolite_id == "compid_46493", 
                            "bilirubin degradation product, C16H18N2O5 (1)*",
                            orig_feature[which(orig_feature$feature_id == metabolite_id),"BIOCHEMICAL"])
  boxplot(dtst[,i] ~ dtst[,"score_group"],
          main = paste( metabolite_name ),
          ylab = expression(paste("Peak intensity (z-scores)")),
          cex.main=0.8 )
}
invisible(dev.off())

####################################################################################
####################################################################################
## identify principle variables
####################################################################################
####################################################################################

pv_res <- iPVs::iPVs(
  metab_full[,metab_list],
  cor_method = "spearman",
  dist_method = "R2",
  hclust_meth = "average",
  cutheight = 0.75
)

pv <- unlist(unname(data.frame(lapply(pv_res$iPV_table$PVs, as.character), stringsAsFactors=FALSE)))

pv_information <- data.frame(matrix(data = pv, nrow = length(pv), ncol = 1, dimnames = list(c(1:length(pv)),"feature_id")))

orig_feature$feature_id <- rownames(orig_feature)
pv_information <- merge(pv_information, orig_feature, by = "feature_id")

pv_id <- pv_information$feature_id
paste0("Number of independent features: ", length(pv_id))

pv_id_with_names <- as.data.frame(pv_id)
pv_id_with_names$feature_name <- orig_feature$BIOCHEMICAL[match(pv_id_with_names$pv_id,orig_feature$feature_id)]
names(pv_id_with_names)[1] <- "feature_id"

# id by biochemical id members of clusters
for (i in 1:length(pv_res$PV_cluster_members)) {
  print(paste0("Cluster ",i))
  cluster_members <- as.data.frame(pv_res$PV_cluster_members[[i]])
  names(cluster_members)[1] <- "feature_id"
  cluster_members$BIOCHEMICAL <- orig_feature$BIOCHEMICAL[match(cluster_members$feature_id,orig_feature$feature_id)]
  print(cluster_members)
}

####################################################################################
####################################################################################
## plot metabolite levels against observed bmi
####################################################################################
####################################################################################

# in this section the relationship between the associated metabolites (from linear regression model) and bmi is evaluated
# first models fitted with both bmi and score group fitted, and then bmi*score group fitted and test for significant interaction done
# finally model run within score group
# all models run on unimputed data with no transformation
# shapiro wilk test done and histograms exported

pdf_out = paste0(data_output_dir,"08.0_histograms_of_rep_metabolite_levels.pdf")
pdf(pdf_out)
par(mfrow=c(4,3))
par(mar = c(2,5,2,2))

result_table <- as.data.frame(matrix(data = NA, nrow = length(metab_list), ncol = 43))
residual_data <- as.data.frame(matrix(data = NA, nrow = dim(dtst)[1], ncol = 1,
                                      dimnames = list(dtst$CLIENT_IDENTIFIER, c("CLIENT_IDENTIFIER"))))
residual_data$CLIENT_IDENTIFIER <- dtst$CLIENT_IDENTIFIER

for (i in 1:length(metab_list)) {
  metabolite_id <- metab_list[i]
  metabolite_name <- orig_feature[which(orig_feature$feature_id == metabolite_id),"BIOCHEMICAL"]
  print(metabolite_id)
  print(metabolite_name)
  # all samples included with bmi predictors (just for power calc)
  fit_temp <- lm(data = dtst,
                 as.formula(paste(metabolite_id,"bmi_F24", sep = "~")))
  print(summary(fit_temp))
  # all samples included with bmi, score group, sex and age as predictors
  fit1 <- lm(data = dtst,
            as.formula(paste(metabolite_id,"bmi_F24 + score_group + sex + age_wks_F24", sep = "~")))
  res <- summary(fit1)
  fit1.robust = lmtest::coeftest(fit1, vcov = sandwich)
  if (i==1) {print(fit1.robust)}
  result_table[i,1] <- metabolite_id
  result_table[i,2] <- metabolite_name
  result_table[i,3:5] <- fit1.robust[2,c(1,2,4)]
  result_table[i,6] <- fit1.robust[2,1] - (1.96*fit1.robust[2,2])
  result_table[i,7] <- fit1.robust[2,1] + (1.96*fit1.robust[2,2])
  result_table[i,8:10] <- fit1.robust[3,c(1,2,4)]
  result_table[i,11:13] <- fit1.robust[4,c(1,2,4)]
  result_table[i,14:16] <- fit1.robust[5,c(1,2,4)]
  result_table[i,17] <- res$adj.r.squared
  result_table[i,18] <- glance(fit1)[5]
  residuals <- as.data.frame(res$residuals, 
                             row.names = rownames(na.omit(dtst[,c("bmi_F24",metabolite_id)])))
  names(residuals) <- metabolite_id
  residuals$CLIENT_IDENTIFIER <- rownames(residuals)
  residual_data <- merge(residual_data, residuals, by = "CLIENT_IDENTIFIER", all.x = TRUE)
  
  # all samples included with bmi, score group, sex and age as predictors plus interaction term between bmi and score group
  fit2 <- lm(data = dtst,
            as.formula(paste(metabolite_id,"bmi_F24 * score_group + sex + age_wks_F24", sep = "~")))
  res <- summary(fit2)
  fit2.robust = lmtest::coeftest(fit2, vcov = sandwich)
  if (i==1) {print(fit2.robust)}
  result_table[i,19:21] <- fit2.robust[2,c(1,2,4)]
  result_table[i,22] <- fit2.robust[2,1] - (1.96*fit2.robust[2,2])
  result_table[i,23] <- fit2.robust[2,1] + (1.96*fit2.robust[2,2])
  result_table[i,24:26] <- fit2.robust[3,c(1,2,4)]
  result_table[i,27:29] <- fit2.robust[6,c(1,2,4)]
  result_table[i,30] <- res$adj.r.squared
  result_table[i,31] <- glance(fit2)[5]
  # check normality
  shapwilk <- shapiro.test(dtst[,metabolite_id])
  result_table[i,32:33] <- c(shapwilk$statistic,shapwilk$p.value) 
  hist(dtst[,metabolite_id],main = metabolite_name,cex.main=0.8)
  # fit model with bmi, score group, sex and age as predictors in low score group
  tmp_dtst_low <- dtst %>% dplyr::filter(score_group == "low")
  fit_low <- lm(data = tmp_dtst_low,
                as.formula(paste(metabolite_id, "bmi_F24 + sex + age_wks_F24", sep = "~")))
  fit_low.robust <- lmtest::coeftest(fit_low, vcov = sandwich)
  result_table[i,34:36] <- fit_low.robust[2,c(1,2,4)]
  if (i==1) {print(fit_low.robust)}
  result_table[i,37] <- fit_low.robust[2,1] - (1.96*fit_low.robust[2,2])
  result_table[i,38] <- fit_low.robust[2,1] + (1.96*fit_low.robust[2,2])
  # fit model with bmi, score group, sex and age as predictors in high score group
  tmp_dtst_high <- dtst %>% dplyr::filter(score_group == "high")
  fit_high <- lm(data = tmp_dtst_high,
                as.formula(paste(metabolite_id, "bmi_F24 + sex + age_wks_F24", sep = "~")))
  fit_high.robust <- lmtest::coeftest(fit_high, vcov = sandwich)
  if (i==1) {print(fit_high.robust)}
  result_table[i,39:41] <- fit_high.robust[2,c(1,2,4)]
  result_table[i,42] <- fit_high.robust[2,1] - (1.96*fit_high.robust[2,2])
  result_table[i,43] <- fit_high.robust[2,1] + (1.96*fit_high.robust[2,2])
}
invisible(dev.off())

names(result_table)[1:18] <- c("feature_id","Metabolite.name","BMI.beta","BMI.Beta.standard.error","BMI.P.value","BMI.Lower.95%.CI","BMI.Upper.95%.CI",
                                "Group.beta","Group.Beta.standard.error","Group.P.value", "Sex.beta", "Sex.Beta.standard.error", "Sex.P.value",
                                "Age.beta","Age.Beta.standard.error", "Age.P.value","Simple.Model.R2", "Simple.Model.P.value")
names(result_table)[19:31] <- c("Interaction.BMI.beta","Interaction.BMI.standard.error","Interaction.BMI.P.value","Interaction.BMI.Lower.95%.CI","Interaction.BMI.Upper.95%.CI",
                                "Interaction.Group.beta","Interaction.Group.Beta.standard.error","Interaction.Group.P.value",
                                "Interaction.BMI*Group.beta","Interaction.BMI*Group.beta.standard.error","Interaction.BMI*Group.P.value","Interaction.Model.R2", "Interaction.Model.P.value")
names(result_table)[32:38] <- c("Shapiro.Wilk.Wstat","shapirowilk_p","LowGroup.BMI.beta", "LowGroup.BMI.beta.standard.error", "LowGroup.BMI.P.value",
                                "LowGroup.BMI.Lower.95%.CI","LowGroup.BMI.Upper.95%.CI")
names(result_table)[39:43] <- c("HighGroup.BMI.beta", "HighGroup.BMI.beta.standard.error", "HighGroup.BMI.P.value",
                                "HighGroup.BMI.Lower.95%.CI","HighGroup.BMI.Upper.95%.CI")

# write results out in full
fwrite(result_table, paste0(data_output_dir, "08.1_metabolite_obs_bmi_regression_results.txt"), sep = "\t")

# write out results for supp table

# Paper Ref: Table S7

# add in missingness stats
metab_miss_percent <- as.data.frame(apply(dtst[,c(2:(ncol(dtst)-6))],2,function(x) (sum(is.na(x))/length(x)))*100 )
metab_miss_percent$feature_id <- rownames(metab_miss_percent)
names(metab_miss_percent)[1] <- "percent_missing"
result_table$`%.missing.values` <- metab_miss_percent$percent_missing[match(result_table$feature_id,metab_miss_percent$feature_id)]

# add flag to indicate which metabolites are representative
result_table$Representative.metabolite <- "No"
result_table[which(result_table$Metabolite.name %in% pv_id_with_names$feature_name), "Representative.metabolite"] <- "Yes"

# add beta from main analysis
result_table$main_results_beta <- merged_res$Beta[match(result_table$feature_id,merged_res$feature_id)]
result_table$Directionally.concordant.with.primary.result <- "No"

for (i in 1:nrow(result_table)){
  main_results_col <- which(names(result_table) == "main_results_beta")
  bmi_beta_col <- which(names(result_table) == "BMI.beta")
  if (result_table[i,main_results_col] >=0 & result_table[i,bmi_beta_col] >= 0) {result_table[i,ncol(result_table)] = "Yes"}
  if (result_table[i,main_results_col] <0 & result_table[i,bmi_beta_col] < 0) {result_table[i,ncol(result_table)] = "Yes"}
}

print("Cross tab of concordance (Yes = 21) and representative metabs (Yes = 13)")
table(result_table$Directionally.concordant.with.primary.result, result_table$Representative.metabolite)

# sort by BMI p
result_table <- result_table[order(result_table$BMI.P.value,decreasing = F),]

# write out version for supplementary tables
fwrite(result_table[,c(2,45,44,32,3:18,47,19:31,34:43)], paste0(data_output_dir, "ForPaper/TableS7_08.1_metabolite_obs_bmi_regression_results.txt"), sep = "\t")

# evaluate shapiro wilk
print("Summarise Shapiro-Wilk stats:")
summary(result_table$Shapiro.Wilk.Wstat)

# plot main beta on bmi beta
pdf_out = paste0(data_output_dir,"08.1_scaled_main_beta_on_bmi_beta.pdf")
pdf(pdf_out)
plot(result_table$BMI.beta,result_table$main_results_beta/2.7)
invisible(dev.off())

####################################################################################
####################################################################################
## plot diagnostics
## check distributions of residuals and their relationship with bmi
####################################################################################
####################################################################################

# Residual QQplots #
pdf_out = paste0(data_output_dir,"08.2_obs_bmi_regression_residuals_qqplot.pdf")
pdf(pdf_out)
par(mfrow=c(4,3))
par(mar = c(2,5,2,2))

for (i in 2:(ncol(residual_data)) ){
  metabolite_id <- names(residual_data)[i]
  metabolite_name <-  orig_feature[which(orig_feature$feature_id == metabolite_id),"BIOCHEMICAL"]
  qqnorm(residual_data[,i],
        main = paste( metabolite_name ),
        cex.main=0.8)
}
invisible(dev.off())

# Residual histograms #
pdf_out = paste0(data_output_dir,"08.3_obs_bmi_regression_residuals_hist.pdf")
pdf(pdf_out)
par(mfrow=c(4,3))
par(mar = c(2,5,2,2))

for (i in 2:(ncol(residual_data)) ){
  metabolite_id <- names(residual_data)[i]
  metabolite_name <- orig_feature[which(orig_feature$feature_id == metabolite_id),"BIOCHEMICAL"]
  hist(residual_data[,i],
         main = paste( metabolite_name ),
         cex.main=0.8)
}
invisible(dev.off())

# Residuals plotted against bmi #
pdf_out = paste0(data_output_dir,"08.4_obs_bmi_regression_residuals_vs_bmi.pdf")
pdf(pdf_out)
par(mfrow=c(4,3))
par(mar = c(2,5,2,2))

for (i in 2:(ncol(residual_data)) ){
  metabolite_id <- names(residual_data)[i]
  metabolite_name <- orig_feature[which(orig_feature$feature_id == metabolite_id),"BIOCHEMICAL"]
  plot(x = residual_data[,i], y = bmi_data$bmi_F24,
       main = paste( metabolite_name ),
       xlab = expression(paste("Residual")),
       ylab = expression(paste("BMI")),
       cex.main=0.8 )
}
invisible(dev.off())


####################################################################################
####################################################################################
## correlation matrix and heatmap for metabolite signals
####################################################################################
####################################################################################

cormat <- round(cor(metab_full[,metab_list],use = "complete.obs"),2)
head(cormat)
melted_cormat <- melt(cormat)
head(melted_cormat)

pdf_out = paste0(data_output_dir,"08.5_assoc_metab_cor_heatmap.pdf")
pdf(pdf_out)
ggplot(data = melted_cormat, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile() +
  theme(axis.text.x=element_text(angle=45,hjust=1),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())
invisible(dev.off())


####################################################################################
####################################################################################
## plot metabolite ~ observed bmi
####################################################################################
####################################################################################

# plot scatters of metabolite on bmi with regression lines
p1 <- plot_pheno_scatter(pv_id, 1, metab_full, pheno_data=bmi_data, orig_feature, result_table, pheno_var = "bmi_F24", pheno_name = "BMI")
p2 <- plot_pheno_scatter(pv_id, 2, metab_full, pheno_data=bmi_data, orig_feature, result_table, pheno_var = "bmi_F24", pheno_name = "BMI")
p3 <- plot_pheno_scatter(pv_id, 3, metab_full, pheno_data=bmi_data, orig_feature, result_table, pheno_var = "bmi_F24", pheno_name = "BMI")
p4 <- plot_pheno_scatter(pv_id, 4, metab_full, pheno_data=bmi_data, orig_feature, result_table, pheno_var = "bmi_F24", pheno_name = "BMI")
p5 <- plot_pheno_scatter(pv_id, 5, metab_full, pheno_data=bmi_data, orig_feature, result_table, pheno_var = "bmi_F24", pheno_name = "BMI")
p6 <- plot_pheno_scatter(pv_id, 6, metab_full, pheno_data=bmi_data, orig_feature, result_table, pheno_var = "bmi_F24", pheno_name = "BMI")
p7 <- plot_pheno_scatter(pv_id, 7, metab_full, pheno_data=bmi_data, orig_feature, result_table, pheno_var = "bmi_F24", pheno_name = "BMI")
p8 <- plot_pheno_scatter(pv_id, 8, metab_full, pheno_data=bmi_data, orig_feature, result_table, pheno_var = "bmi_F24", pheno_name = "BMI")
p9 <- plot_pheno_scatter(pv_id, 9, metab_full, pheno_data=bmi_data, orig_feature, result_table, pheno_var = "bmi_F24", pheno_name = "BMI")
p10 <- plot_pheno_scatter(pv_id, 10, metab_full, pheno_data=bmi_data, orig_feature, result_table, pheno_var = "bmi_F24", pheno_name = "BMI")
p11 <- plot_pheno_scatter(pv_id, 11, metab_full, pheno_data=bmi_data, orig_feature, result_table, pheno_var = "bmi_F24", pheno_name = "BMI")
p12 <- plot_pheno_scatter(pv_id, 12, metab_full, pheno_data=bmi_data, orig_feature, result_table, pheno_var = "bmi_F24", pheno_name = "BMI")
p13 <- plot_pheno_scatter(pv_id, 13, metab_full, pheno_data=bmi_data, orig_feature, result_table, pheno_var = "bmi_F24", pheno_name = "BMI")
p14 <- plot_pheno_scatter(pv_id, 14, metab_full, pheno_data=bmi_data, orig_feature, result_table, pheno_var = "bmi_F24", pheno_name = "BMI")
p15 <- plot_pheno_scatter(pv_id, 15, metab_full, pheno_data=bmi_data, orig_feature, result_table, pheno_var = "bmi_F24", pheno_name = "BMI")


# Paper Ref: Fig 5

p_patchwork <- (p15|p2|p6)/(p14|p13|p1|p9)
p_patchwork
ggsave(filename = paste0(data_output_dir,"ForPaper/Fig5_08.5a_obs_bmi_scatters_main.pdf"),
       plot = last_plot(),width = 40,height = 20,units = c("cm"),dpi = 320)

# Paper Ref: Fig S4

p_patchwork <- (p3|p4|p5|p7)/(p8|p10|p11|p12)
p_patchwork
ggsave(filename = paste0(data_output_dir,"ForPaper/FigS5_08.5b_obs_bmi_scatters_supp.pdf"),
       plot = last_plot(),width = 35,height = 20,units = c("cm"),dpi = 320)


####################################################################################
####################################################################################
## Compare primary results to BMI observational
####################################################################################
####################################################################################

pdf_out = paste0(data_output_dir,"ForPaper/FigS4_08.6_obs_beta_on_primary_beta.pdf")
pdf(file=pdf_out)
plot(result_table$main_results_beta,result_table$BMI.beta, xlab = "Model 1 (BMI GRS group) beta", ylab = "Observational BMI beta", pch = 4)
abline(h=0,v=0,lty=3)
invisible(dev.off())

#################
###### END ######
#################

print(sessionInfo())

q()
