# script to prepare phenotype (not metabolite) data


####################
###### SET UP ######
####################

# read in parameter file (specified on command line)
source("parameter_files/parameters_for_r.R")

# move to working directory 
setwd(working_dir)

# load library
library(dplyr)

#####################
###### DATA IN ######
#####################

# read in metabolite sample metadata
metab_sample <- read.table(file=paste0(metab_input_dir,"B3194_BMI_RbG_2020_02_03_Metabolon_sampledata.txt"), h=T, stringsAsFactors=F)
print(paste0("There are ", nrow(metab_sample), " samples in the Metabolon sample data file."))

# read in pheno data used to do sample selection (RbG)
pheno_data <- read.table(file="input/selection/sampled.txt",sep="\t",h=T)
print(paste0("There are ", nrow(pheno_data), " individuals in the selection input file for RbG (from age 24 clinic)."))

# read in pheno data used to do sample selection (RbBMI)
pheno_data_bmi <- read.csv(file = paste0(data_intermediate_dir,"bmi_rbg_meta_processed_with_imputed_weight.csv"))
pheno_data_bmi <- pheno_data_bmi[,c("aln","qlet","bmi_F24")]
print(paste0("There are ", nrow(pheno_data_bmi), " individuals in the selection input file for RbBMI."))

##########################
###### PROCESS DATA ######
########## RbG ###########
##########################

### restrict pheno data to those selected for study
# replace NA with 0's (so the "included" var can be used for subsetting)
pheno_data[["included"]][is.na(pheno_data[["included"]])] <- 0
# restrict to selected
selected <- pheno_data[pheno_data$included == 1,]
print(paste0("There are ", nrow(selected), " individuals in the selected subset."))

### merge 
# make common identifier for merging the metabolite sample data with selection info.
selected$aln_qlet <- paste(selected$aln,selected$qlet,sep="_")

# identify high/low BMI PRS groups
mean_yengo <- mean(selected$yengo_score)
selected$group <- "NA"
# define high and low groups based on whether an individual's score was above or below the mean score
selected[which(selected$yengo_score >= mean_yengo), "group"] <- "high"
selected[which(selected$yengo_score < mean_yengo), "group"] <- "low"

# import group to metab meta file
orig_sample <- metab_sample
orig_sample$score_group <- selected$group[match(orig_sample$CLIENT_IDENTIFIER,selected$aln_qlet)]

# identify duplicates (n=4 pairs)
n_occur <- data.frame(table(orig_sample$CLIENT_IDENTIFIER))
print(paste0("The following individuals had duplicate samples (two) analysed:"))
n_occur[n_occur$Freq > 1,]

# identify those selected but without metabolite data (n=4)
missing_ids <- which(!(selected$aln_qlet %in% orig_sample$CLIENT_IDENTIFIER))
print(paste0("The following individuals had no samples in analysed:"))
selected[missing_ids,]

######################
###### SAVE OUT ######
######################

# save out
write.table(orig_sample,file = paste0(data_intermediate_dir,"Metabolon_sampledata_with_group.txt"), sep="\t", row.names=F, quote=F)

##########################
###### PROCESS DATA ######
######### RbBMI ##########
##########################

# calculate BMI range
# check class is numeric
class(pheno_data_bmi$bmi_F24)
pheno_data_bmi$bmi_F24[which(pheno_data_bmi$bmi_F24 < 0)] <- NA
bmi_range <- quantile(pheno_data_bmi$bmi_F24, probs = c(0.3,0.7), na.rm = T)

# select individuals within BMI range of interest
bmi_low <- pheno_data_bmi %>% filter(bmi_F24 < bmi_range[1])
bmi_high <- pheno_data_bmi %>% filter(bmi_F24 > bmi_range[2])

# randomly select 380 individuals in each group
set.seed(1234)
bmi_low <- bmi_low[sample(nrow(bmi_low),380, replace=F),]
bmi_high <- bmi_high[sample(nrow(bmi_high),380, replace=F),]

# add group labels
bmi_low$bmi_group <- "low"
bmi_high$bmi_group <- "high"

selected_bmi <- rbind(bmi_low,bmi_high)
selected_bmi$aln_qlet <- paste(selected_bmi$aln,selected_bmi$qlet,sep="_")

######################
###### SAVE OUT ######
######################

# save out
write.table(selected_bmi,file = paste0(data_intermediate_dir,"BMI_sampledata_with_group.txt"), 
            sep="\t", row.names=F, quote=F)

#################
###### END ######
#################

print(sessionInfo())

q()
