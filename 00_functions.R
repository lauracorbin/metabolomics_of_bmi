# functions  for BMI RbG analysis
# functions script sourced by other scripts in project

## load libraries
library(directlabels)

# winsorise function
winsorize_x = function(x, cut=0.01) {
  #winsorize_x = function(x, nsd=5) {
  cut_point_top <- quantile(x,1-cut,na.rm=T)
  cut_point_bottom <- quantile(x,cut,na.rm=T)
  #cut_point_top <- mean(x,na.rm=T) + (nsd*(sd(x,na.rm=T)))
  #cut_point_bottom <- mean(x,na.rm=T) - (nsd*(sd(x,na.rm=T)))
  i = which(x >= cut_point_top)
  x[i] = cut_point_top
  j = which(x <= cut_point_bottom)
  x[j] = cut_point_bottom
  return(x)
}


# convert to pa
pa.convert = function(dtst=dtst) {
  output <- apply(dtst,2,function(x) ifelse(is.na(x),0,1))
  metab_number <- (ncol(output)-10)/2
  for (j in 1:nrow(output)) {
    # set to NA if entire baseline sample missing
    if ( sum(output[j,1:metab_number] == 0) == metab_number ) {
      output[j,1:metab_number] <- NA
    }
  }
  return(as.data.frame(output))
}

# linear analysis function (simple)
linear_analysis_simple <- function(data, var_list, covar_metab) {
  dtst <- merge(data,covar_metab,by.x = "CLIENT_IDENTIFIER", by.y = "aln_qlet")   ## define dataset 
  linear_results = as.data.frame(matrix(data = NA, nrow = length(var_list), ncol = 12))
  if (levels(dtst$score_group)[1] == "high") {
    dtst <- within(dtst,score_group <- relevel(score_group, ref = "low"))
  }
  s = base::seq(from=1,to=length(feature_list_low_miss), by=100)
  # modelling loop
  for (i in 1:length(var_list)){
    if(i %in% s){
      print(paste0( "Now processing metabolite ", i, " of ", length(var_list) ))
    }
    mtb_name <- var_list[i]
    linear_results[i,1] <- mtb_name
    mtb_col_ref <- which(colnames(dtst) == mtb_name)
    # subset data to single feature
    mtb <- dtst[!is.na(dtst[mtb_col_ref]),c("SAMPLE_NAME","score_group",mtb_name)]
    # record no. of observations
    linear_results[i,2] <- nrow(mtb)
    # check numbers (important if using unimputed data)
    group_fac <- length(unique(mtb$score_group))
    if ( (group_fac > 1) & (nrow(mtb) > 10) ) {
      # rank transform
      mtb[,3] <- rntransform(mtb[,3])
      linear_results[i,3] <- "simple"
      fit <- lm(as.formula(paste0(mtb_name,"~ score_group")), 
                data=na.omit(mtb))
      r2 <- summary(fit)$r.squared
      linear_results[i,4] <- r2
      coef = summary(fit)$coefficients
      group_coef = coef[2 ,]
      a = Anova(fit)
      eta = a[,1]/sum(a[,1])
      linear_results[i,5:8] <- group_coef
      linear_results[i,9] <- group_coef[1] - (group_coef[2]*1.96)
      linear_results[i,10] <- group_coef[1] + (group_coef[2]*1.96)
      linear_results[i,11:12] <- eta
    } else {
      linear_results[i,5:12] <- NA
    }
  }
  names(linear_results)[1] <- "feature_id"
  names(linear_results)[2] <- "n_samples"
  names(linear_results)[3] <- "model_type"
  names(linear_results)[4] <- "model_total_r2"
  names(linear_results)[5] <- "lm_group_beta"
  names(linear_results)[6] <- "lm_group_SE"
  names(linear_results)[7] <- "lm_group_tval"
  names(linear_results)[8] <- "lm_group_p"
  names(linear_results)[9] <- "lm_group_lowerCI"
  names(linear_results)[10] <- "lm_group_upperCI"
  names(linear_results)[11:12] = paste0("eta2_", rownames(a))
  
  ## adjust p-val
  linear_results$lm_group_adjp = p.adjust(linear_results[, "lm_group_p"], method = "BH")
  
  ## identify significant features
  linear_results$linear_flag = 0
  linear_results[which(linear_results$lm_group_adjp <= 0.05),"linear_flag"] <- 1

  # sort by p
  linear_results <- linear_results[order(linear_results$lm_group_p,decreasing = F),]
  
  # add rank
  linear_results$linear_score_p_rank <- rank(linear_results$lm_group_p)
  
  # save out result
  return(linear_results)
}

# linear analysis function (complex)
linear_analysis_complex <- function(data, var_list, covar_metab) {
  dtst <- merge(data,covar_metab,by.x = "CLIENT_IDENTIFIER", by.y = "aln_qlet")   ## define dataset 
  linear_results = as.data.frame(matrix(data = NA, nrow = length(var_list), ncol = 14))
  if (levels(dtst$score_group)[1] == "high") {
    dtst <- within(dtst,score_group <- relevel(score_group, ref = "low"))
  }
  s = base::seq(from=1,to=length(feature_list_low_miss), by=100)
  # modelling loop
  for (i in 1:length(var_list)){
    if(i %in% s){
      print(paste0( "Now processing metabolite ", i, " of ", length(var_list) ))
    }
    mtb_name <- var_list[i]
    linear_results[i,1] <- mtb_name
    mtb_col_ref <- which(colnames(dtst) == mtb_name)
    # subset data to single feature
    mtb <- dtst[!is.na(dtst[mtb_col_ref]),c("SAMPLE_NAME","score_group",mtb_name,"mat_social_class","pat_social_class")]
    # record no. of observations
    linear_results[i,2] <- nrow(mtb)
    # check numbers
    group_fac <- length(unique(mtb$score_group))
    if ( (group_fac > 1) & (nrow(mtb) > 10) ) {
      # rank transform
      mtb[,3] <- rntransform(mtb[,3])
      linear_results[i,3] <- "complex"
      fit <- lm(as.formula(paste0(mtb_name,"~ score_group + pat_social_class + mat_social_class")), 
                data=na.omit(mtb))
      r2 <- summary(fit)$r.squared
      linear_results[i,4] <- r2
      coef = summary(fit)$coefficients
      group_coef = coef[2 ,]
      a = Anova(fit)
      eta = a[,1]/sum(a[,1])
      linear_results[i,5:8] <- group_coef
      linear_results[i,9] <- group_coef[1] - (group_coef[2]*1.96)
      linear_results[i,10] <- group_coef[1] + (group_coef[2]*1.96)
      linear_results[i,11:14] <- eta
    } else {
      linear_results[i,5:14] <- NA
    }
  }
  names(linear_results)[1] <- "feature_id"
  names(linear_results)[2] <- "n_samples"
  names(linear_results)[3] <- "model_type"
  names(linear_results)[4] <- "model_total_r2"
  names(linear_results)[5] <- "lm_group_beta"
  names(linear_results)[6] <- "lm_group_SE"
  names(linear_results)[7] <- "lm_group_tval"
  names(linear_results)[8] <- "lm_group_p"
  names(linear_results)[9] <- "lm_group_lowerCI"
  names(linear_results)[10] <- "lm_group_upperCI"
  names(linear_results)[11:14] = paste0("eta2_", rownames(a))
  
  ## adjust p-val
  linear_results$lm_group_adjp = p.adjust(linear_results[, "lm_group_p"], method = "BH")
  
  ## identify significant features
  linear_results$linear_flag = 0
  linear_results[which(linear_results$lm_group_adjp <= 0.05),"linear_flag"] <- 1

  # sort by p
  linear_results <- linear_results[order(linear_results$lm_group_p,decreasing = F),]
  
  # add rank
  linear_results$linear_score_p_rank <- rank(linear_results$lm_group_p)
  
  # save out result
  return(linear_results)
}

# linear analysis function (BMI model)
linear_analysis_BMI <- function(data, var_list, covar_metab) {
  dtst <- merge(data,covar_metab,by.x = "CLIENT_IDENTIFIER", by.y = "aln_qlet")   ## define dataset 
  linear_results = as.data.frame(matrix(data = NA, nrow = length(var_list), ncol = 17))
  if (levels(dtst$score_group)[1] == "high") {
    dtst <- within(dtst,score_group <- relevel(score_group, ref = "low"))
  }
  s = base::seq(from=1,to=length(feature_list_low_miss), by=100)
  # modelling loop
  for (i in 1:length(var_list)){
    if(i %in% s){
      print(paste0( "Now processing metabolite ", i, " of ", length(var_list) ))
    }
    mtb_name <- var_list[i]
    linear_results[i,1] <- mtb_name
    mtb_col_ref <- which(colnames(dtst) == mtb_name)
    # subset data to single feature
    mtb <- dtst[!is.na(dtst[mtb_col_ref]),c("SAMPLE_NAME","score_group",mtb_name,
                                            "mat_social_class","pat_social_class","bmi_F24","sex","age_wks_F24")]
    # record no. of observations
    linear_results[i,2] <- nrow(mtb)
    # check numbers
    group_fac <- length(unique(mtb$score_group))
    if ( (group_fac > 1) & (nrow(mtb) > 10) ) {
      # rank transform
      mtb[,3] <- rntransform(mtb[,3])
      linear_results[i,3] <- "BMI"
      fit <- lm(as.formula(paste0(mtb_name,"~ bmi_F24 + sex + age_wks_F24 + score_group + pat_social_class + mat_social_class")), 
                data=na.omit(mtb))
      r2 <- summary(fit)$r.squared
      linear_results[i,4] <- r2
      coef = summary(fit)$coefficients
      group_coef = coef[2 ,]
      a = Anova(fit)
      eta = a[,1]/sum(a[,1])
      linear_results[i,5:8] <- group_coef
      linear_results[i,9] <- group_coef[1] - (group_coef[2]*1.96)
      linear_results[i,10] <- group_coef[1] + (group_coef[2]*1.96)
      linear_results[i,11:17] <- eta
    } else {
      linear_results[i,5:17] <- NA
    }
  }
  names(linear_results)[1] <- "feature_id"
  names(linear_results)[2] <- "n_samples"
  names(linear_results)[3] <- "model_type"
  names(linear_results)[4] <- "model_total_r2"
  names(linear_results)[5] <- "lm_BMI_beta"
  names(linear_results)[6] <- "lm_BMI_SE"
  names(linear_results)[7] <- "lm_BMI_tval"
  names(linear_results)[8] <- "lm_BMI_p"
  names(linear_results)[9] <- "lm_BMI_lowerCI"
  names(linear_results)[10] <- "lm_BMI_upperCI"
  names(linear_results)[11:17] = paste0("eta2_", rownames(a))
  
  ## adjust p-val
  linear_results$lm_BMI_adjp = p.adjust(linear_results[, "lm_BMI_p"], method = "BH")
  
  ## identify significant features
  linear_results$linear_flag = 0
  linear_results[which(linear_results$lm_BMI_adjp <= 0.05),"linear_flag"] <- 1

  # sort by p
  linear_results <- linear_results[order(linear_results$lm_BMI_p,decreasing = F),]
  
  # add rank
  linear_results$linear_score_p_rank <- rank(linear_results$lm_BMI_p)
  
  # save out result
  return(linear_results)
}

# Q-Q plot
gg_qqplot <- function(ps, feature_id, ci = 0.95) {
  z = qnorm(ps / 2)
  lambda = round(median(z^2) / qchisq(0.5, 1), 3)
  n  <- length(ps)
  df <- data.frame(
    id = feature_id,
    observed = -log10(sort(ps)),
    expected = -log10(ppoints(n)),
    clower   = -log10(qbeta(p = (1 - ci) / 2, shape1 = 1:n, shape2 = n:1)),
    cupper   = -log10(qbeta(p = (1 + ci) / 2, shape1 = 1:n, shape2 = n:1))
  )
  log10Pe <- expression(paste("Expected -log"[10], plain(P)))
  log10Po <- expression(paste("Observed -log"[10], plain(P)))
  ggplot(df) +
    geom_point(aes(expected, observed), shape = 1, size = 3) +
    geom_abline(intercept = 0, slope = 1, alpha = 0.5) +
    geom_line(aes(expected, cupper), linetype = 2) +
    geom_line(aes(expected, clower), linetype = 2) +
    xlab(log10Pe) +
    ylab(log10Po) +
    labs(subtitle = paste0("lambda = ",lambda)) +
    theme_bw() +
    geom_point(data=df[which(df$id %in% c("Glucose_F24", "Glc_F24")),],
               aes(expected, observed), colour="red", size=3)
}

# histogram and QQ plot of pval of linear analysis
linear_analysis_plots <- function(data, model) {
  p_hist <- ggplot(data = data, aes(x=lm_group_p)) + 
    geom_histogram(bins = 30) + 
    labs(title = paste0("Distribution of unadjusted p-value (",model," model)"),
         x = "Unadjusted p-value",y = "Frequency") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5))
  p_qq <- gg_qqplot(data$lm_group_p, data$feature_id)
  
  return(list(p_hist,p_qq))
}

linear_analysis_plots_BMI <- function(data, model) {
  p_hist <- ggplot(data = data, aes(x=lm_BMI_p)) + 
    geom_histogram(bins = 30) + 
    labs(title = paste0("Distribution of unadjusted p-value (",model," model)"),
         x = "Unadjusted p-value",y = "Frequency") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5))
  p_qq <- gg_qqplot(data$lm_BMI_p, data$feature_id)
  
  return(list(p_hist,p_qq))
}


# make metabolite adiposity scatters plots
plot_pheno_scatter <- function(pv_id, n, metab_full, pheno_data, orig_feature, result_table, pheno_var, pheno_name) {
  metabolite_id <- as.character(pv_id[n])
  metabolite_name <- ifelse(metabolite_id == "compid_46493", 
                            "bilirubin degradation product, C16H18N2O5 (1)",
                            orig_feature[which(orig_feature$feature_id == metabolite_id),"BIOCHEMICAL"])
  plot.dtst <- merge(metab_full[,c("CLIENT_IDENTIFIER",metabolite_id)], pheno_data, 
                by.x = "CLIENT_IDENTIFIER", by.y = "aln_qlet")
  plot.dtst <- na.omit(plot.dtst)
  colnames(plot.dtst)[2] <- "metabolite_level"
  plot.dtst$metabolite_level <- scale(plot.dtst$metabolite_level)
  
  yrng <- range(plot.dtst$metabolite_level)
  xrng <- range(plot.dtst[,`pheno_var`])
  beta <- round(result_table[which(result_table$feature_id == metabolite_id),paste0(pheno_name,".beta")], digits = 3)
  pval <- round(result_table[which(result_table$feature_id == metabolite_id),paste0(pheno_name,".P.value")], digits = 3)
  pval_int <- round(result_table[which(result_table$feature_id == metabolite_id),paste0("Interaction.",pheno_name,"*Group.P.value")], digits = 3)
  lci <- round(result_table[which(result_table$feature_id == metabolite_id),paste0(pheno_name,".Lower.95%.CI")], digits = 3)
  uci <- round(result_table[which(result_table$feature_id == metabolite_id),paste0(pheno_name,".Upper.95%.CI")], digits = 3)
  
  # define x axis label
  if (pheno_name == "BMI") {xlabel <- expression(paste("Measured BMI (kg/",m ^ 2,")"))}
  if (pheno_name == "lean_mass") {xlabel <- "Total body lean mass (kg)"}
  if (pheno_name == "lean_ratio") {xlabel <- "Total body lean mass : Total fat mass ratio"}
  
  if (pval >= 0.05 | pval_int < 0.05) {
    
    beta_high <- round(result_table[which(result_table$feature_id == metabolite_id),paste0("HighGroup.",pheno_name,".beta")], digits = 3)
    pval_high <- round(result_table[which(result_table$feature_id == metabolite_id),paste0("HighGroup.",pheno_name,".P.value")], digits = 3)
    beta_low <- round(result_table[which(result_table$feature_id == metabolite_id),paste0("LowGroup.",pheno_name,".beta")], digits = 3)
    pval_low <- round(result_table[which(result_table$feature_id == metabolite_id),paste0("LowGroup.",pheno_name,".P.value")], digits = 3)
    lci_high <- round(result_table[which(result_table$feature_id == metabolite_id),paste0("HighGroup.",pheno_name,".Lower.95%.CI")], digits = 3)
    uci_high <- round(result_table[which(result_table$feature_id == metabolite_id),paste0("HighGroup.",pheno_name,".Upper.95%.CI")], digits = 3)
    lci_low <- round(result_table[which(result_table$feature_id == metabolite_id),paste0("LowGroup.",pheno_name,".Lower.95%.CI")], digits = 3)
    uci_low <- round(result_table[which(result_table$feature_id == metabolite_id),paste0("LowGroup.",pheno_name,".Upper.95%.CI")], digits = 3)
    
    p <- ggplot(plot.dtst, aes(x = plot.dtst[,`pheno_var`], 
                          y = metabolite_level, 
                          color = score_group)) +
      geom_point(alpha = 0.3) +
      stat_smooth(aes(fill = score_group, color = score_group), method = "lm") +
      theme_bw() +
      theme(legend.position = "none",
            axis.title.x=element_text(size=8),
            axis.title.y=element_text(size=8),
            axis.text.y=element_text(hjust=0.5),
            plot.subtitle = element_text(size = 8)) +
      xlab(xlabel) +
      ylab("Scaled metabolite level") +
      labs(subtitle = paste0(metabolite_name)) +
      annotate(geom = "text" ,size=3, x = ( xrng[2]-xrng[1] )*0.5 + xrng[1], y = yrng[2],
               label = paste("list(~beta[high.BMI.GRS]~(CI[95~'%']) == ", beta_high,"~(" ,"list(", lci_high, "," , uci_high,  ")", ")",")", sep = ""), parse = TRUE, vjust = 0.5) +
      annotate(geom = "text" ,size=3, x = ( xrng[2]-xrng[1] )*0.5 + xrng[1], y = (yrng[2]-yrng[1])*0.93 + yrng[1],
             label = paste("list(~beta[low.BMI.GRS]~(CI[95~'%']) == ", beta_low,"~(" ,"list(", lci_low, "," , uci_low,  ")", ")" , ")", sep = "")  , parse = TRUE, vjust = 0.5) +
      scale_color_jco() + scale_fill_jco()
  } else {
    p <- ggplot(plot.dtst, aes(x = plot.dtst[,`pheno_var`], 
                          y = metabolite_level, 
                          color = score_group)) +
      geom_point(alpha = 0.3) +
      stat_smooth(aes(fill = score_group, color = score_group), method = "lm") +
      theme_bw() +
      theme(legend.position = "none",
            axis.title.x=element_text(size=8),
            axis.title.y=element_text(size=8),
            axis.text.y=element_text(hjust=0.5),
            plot.subtitle = element_text(size = 8)) +
      xlab(xlabel) +
      ylab("Scaled metabolite level") +
      labs(subtitle = paste0(metabolite_name)) +
      annotate(geom = "text" ,size=3, x = ( xrng[2]-xrng[1] )*0.5 + xrng[1], y = yrng[2],
               label = paste("list(~beta[overall]~(CI[95~'%']) == ", beta,"~(" ,"list(", lci, "," , uci,  ")", ")" ,")", sep = "")  , parse = TRUE, vjust = 0.5) +
      scale_color_jco() + scale_fill_jco()
  }
  return(p)

}

# function for scientific notation (for expression)
# adapted from: 
# https://stackoverflow.com/questions/11610377/how-do-i-change-the-formatting-of-numbers-on-an-axis-with-ggplot
fancy_scientific <- function(l) {
  # turn in to character string in scientific notation
  l <- format(l,digits = 3, scientific = TRUE)
  # quote the part before the exponent to keep all the digits
  l <- gsub("^(.*)e", "'\\1'e", l)
  # turn the 'e+' into plotmath format
  l <- gsub("e", "%*%10^", l)
  return(l)
}



