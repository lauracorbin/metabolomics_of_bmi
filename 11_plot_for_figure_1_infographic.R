# script for making infographic plots

####################
###### SET UP ######
####################

# read in parameter file (specified on command line)
source("parameter_files/parameters_for_r.R")

# move to working directory 
setwd(working_dir)

## load libraries
library(dplyr)
library(ggplot2)

####################################################################################
####################################################################################
## read in data
####################################################################################
####################################################################################

phenofile_full <- read.table(file = paste0(data_intermediate_dir,"bmi_rbg_meta_processed_with_imputed_weight.csv"), h=T, sep=",", stringsAsFactors = F)
sample_data <- read.table(file="input/selection/sampled.txt",sep="\t",h=T)
score_data <- read.table("../selection01/Wade_20200317/yengo_score.profile",header=T)

####################################################################################
####################################################################################
## make plot 2 - adulthood
####################################################################################
####################################################################################

data <- merge(sample_data, phenofile_full[,c("aln","qlet","bmi_F24")], by = c("aln","qlet"))
mean_score <- mean(data$yengo_score, na.rm = T)
data <- data %>% dplyr::mutate(group = ifelse(!is.na(included), ifelse(yengo_score > mean_score, "high", "low"), "others"))
data <- data %>% dplyr::filter(bmi_F24 > 0)
data$aln_qlet <- paste0(data$aln, data$qlet)
data <- data[,c(-1,-2)]

data_full <- merge(score_data, data, by.x="FID", by.y="aln_qlet")

cols <- c("low" = "#EFC000FF", "high" = "#0073C2FF", "others" = "#868686FF")

ggplot(data=data_full, aes(x = SCORESUM, y = bmi_F24)) +
  geom_point(aes(color = factor(group)), alpha = 0.5, size = 1) +
  labs(x = "BMI Genetic Risk Score",
       y = "BMI at 24 Years") +
  geom_smooth(method='lm', color = "#003C67FF") +
  scale_color_manual(values=cols,
                     breaks = c("low", "high")) +
  theme_classic() +
  theme(legend.title = element_blank(),
        legend.justification=c(0,1), 
        legend.position=c(0.05, 1.0),
        legend.background = element_blank(),
        legend.key = element_blank()) +
  guides(color = guide_legend(override.aes = list(size = 3)))

ggsave(paste0(data_output_dir,"ForPaper/Fig1B_adulthood.png"), 
       dpi = 320, width=3.5, height=2)

####################################################################################
####################################################################################
## make plot 1 - conception
####################################################################################
####################################################################################

ggplot(data=data_full, aes(x = SCORESUM, fill = factor(group))) +
  geom_histogram(binwidth = 0.1, position="stack", alpha = 0.7) +
  scale_fill_manual(values=cols,
                     breaks = c("low", "high")) +
  labs(x = "BMI Genetic Risk Score",
       y = "Frequency") +
  theme_classic() +
  theme(legend.title = element_blank(),
        legend.justification=c(0,1), 
        legend.position=c(0.05, 1.0),
        legend.background = element_blank(),
        legend.key = element_blank()) +
  guides(color = guide_legend(override.aes = list(size = 3)))

ggsave(paste0(data_output_dir,"ForPaper/Fig1A_conception.png"),
       dpi = 320, width=3.5, height=2)


#################
###### END ######
#################

print(sessionInfo())

q()







