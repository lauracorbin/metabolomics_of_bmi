# script for sample quality control
# generate a list of duplicated samples to remove
# and a list of individuals to remove in confounder analysis

####################
###### SET UP ######
####################

# read in parameter file (specified on command line)
source("parameter_files/parameters_for_r.R")

# move to working directory 
setwd(working_dir)

# load library
library(dplyr)

####################################################################################
####################################################################################
## read in data
####################################################################################
####################################################################################

## Metabolon data ##
## original scale version (after QC)
orig_data <- read.table(paste0(metab_input_dir,"qc_data/B3194_BMI_RbG_2020_02_03_QCd_data.txt"), h=T)
## sample metadata with group added
orig_sample <- read.table(paste0(data_intermediate_dir,"Metabolon_sampledata_with_group.txt"), 
                          h=T, stringsAsFactors = F)

## read in list of ALSPAC IDs to remove because they're part of a quad set
quad_exc_list <- read.table(paste0(data_intermediate_dir, "quad_exc_list.txt"), h=F, stringsAsFactors = F)

## read in list of ALSPAC IDs to remove because they have subsequently withdrawn consent
consent_withdrawn_list <- read.table(paste0(data_intermediate_dir, "consent_withdrawn_list.txt"), h=F, stringsAsFactors = F)

####################################################################################
####################################################################################
## generate list of samples to remove
####################################################################################
####################################################################################

# bring in sample id
orig_data$SAMPLE_NAME <- rownames(orig_data)

# check duplicates
duplicated_ids <- orig_sample[which(duplicated(orig_sample$CLIENT_IDENTIFIER)),]$CLIENT_IDENTIFIER
duplicated_ids

# find duplicated sample with higher missingness
metab_dups <- orig_sample %>% filter(CLIENT_IDENTIFIER %in% duplicated_ids)
metab_dups_full <- merge(metab_dups, orig_data, by = "SAMPLE_NAME")
metab_dups$sample_missingness <- apply(metab_dups_full[,names(orig_data)],1,function(x) sum(is.na(x))/length(x))

samples_to_remove <- c()
for (i in 1:length(duplicated_ids)) {
  id <- duplicated_ids[i]
  tmp <- metab_dups %>% filter(CLIENT_IDENTIFIER %in% id)
  samples_to_remove <- c(samples_to_remove, tmp[order(tmp$sample_missingness),]$SAMPLE_NAME[-1])
}

# save out
saveRDS(samples_to_remove, file = paste0(data_intermediate_dir,"samples_to_remove.rds"))

####################################################################################
####################################################################################
## generate list of individuals to remove
####################################################################################
####################################################################################

# define "not in" operator
`%notin%` <- Negate(`%in%`)

# find individuals removed during preliminary QC by MetaboQC
# compare individuals in the QC'd metbolite data file to those in the original sample list
removed_samples <- orig_sample %>% filter(SAMPLE_NAME %notin% orig_data$SAMPLE_NAME)
individuals_to_remove <- removed_samples$CLIENT_IDENTIFIER
# add those individuals that are part of a quad group to the exclusion list
individuals_to_remove <- c(individuals_to_remove, quad_exc_list$V1, consent_withdrawn_list$V1)

saveRDS(individuals_to_remove, file = paste0(data_intermediate_dir,"individuals_to_remove.rds"))


#################
###### END ######
#################

print(sessionInfo())

q()
