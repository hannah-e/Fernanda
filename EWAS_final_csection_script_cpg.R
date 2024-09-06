#######################################################################-
#        EWAS DATA PREP: C-SECTION - OFFSPRING DNA METHYLATION         #-
#######################################################################-
#last modified: 20 August 2024

#---------------------------------------------------------------------#
#                             HOUSEKEEPING                            ----
#---------------------------------------------------------------------#

rm(list=ls()) #Remove any existing objects in R 

#set working directory
setwd("/user/work/dy20206/ALSPAC/C_Section")

#load packages
if (!require("pacman")) install.packages("pacman")
pacman::p_load("haven", "dplyr", "purrr", "R.utils", "data.table", "plyr", "matrixStats", "MASS")
pacman::p_load_gh("perishky/meffil", "perishky/ewaff")

#check packages are loaded
(.packages())

# Set results directory
workdir <- "/user/work/dy20206/ALSPAC/C_Section/results_extra"

# Set the cohort info
cohort <- "ALSPAC" # name of your cohort
date <- format(Sys.time(), "%Y-%m-%d")
analyst <- "FMB" # initials of the analyst

#---------------------------------------------------------------------#
#                             LOAD DATA                               ----
#---------------------------------------------------------------------#

#load beta matrix
load("Dataset_final_PACE_csection_cord_F7_v2.rda") #v2 removes any participants missing data on csection (so results here are consistent with those presented for m10a-f)

#---------------------------------------------------------------------#
#      CREATE SURROGATE VARIABLES TO ACCOUNT FOR TECHNICAL BATCH      ----
#---------------------------------------------------------------------#
#we need to do this for each model separately, ethnicity isn't included because all ethnic groups = European

#create lists of covariates
#main models:
  #at birth
covariates_m1 <- c("smoking_mum",
                   "male_sex",
                   "education_mum",
                   "nulliparous",
                   "gestational_age",
                   "birthweight",
                   "bmi_mum") 
covariates_m2 <- c("smoking_mum",
                   "male_sex",
                   "education_mum",
                   "nulliparous",
                   "gestational_age",
                   "birthweight",
                   "bmi_mum",
                   "Bcell_cord", "CD14_cord", "CD4T_cord", "CD8T_cord", "Gran_cord", "NK_cord")
  #childhood
covariates_m3 <- c("smoking_mum",
                   "male_sex",
                   "education_mum",
                   "nulliparous",
                   "age_child",
                   "birthweight",
                   "bmi_mum") 
covariates_m4 <- c("smoking_mum",
                   "male_sex",
                   "education_mum",
                   "nulliparous",
                   "age_child",
                   "birthweight",
                   "bmi_mum",
                   "Bcell_child", "Mono_child", "CD4T_child", "CD8T_child", "Gran_child", "NK_child") 
#sensitivity models:
  #at birth
covariates_m5 <- c("smoking_mum",
                   "male_sex",
                   "education_mum",
                   "nulliparous",
                   "gestational_age",
                   "birthweight",
                   "bmi_mum")
covariates_m6 <- c("smoking_mum",
                   "male_sex",
                   "education_mum",
                   "nulliparous",
                   "gestational_age",
                   "birthweight",
                   "bmi_mum",
                   "Bcell_cord", "CD14_cord", "CD4T_cord", "CD8T_cord", "Gran_cord", "NK_cord")
  #childhood
covariates_m7 <- c("smoking_mum",
                   "male_sex",
                   "education_mum",
                   "nulliparous",
                   "age_child",
                   "birthweight",
                   "bmi_mum")
covariates_m8 <- c("smoking_mum",
                   "male_sex",
                   "education_mum",
                   "nulliparous",
                   "age_child",
                   "birthweight",
                   "bmi_mum",
                   "Bcell_child", "Mono_child", "CD4T_child", "CD8T_child", "Gran_child", "NK_child")
  #cell type models
  #at birth
covariates_m9 <- c("smoking_mum",
                   "male_sex",
                   "education_mum",
                   "nulliparous",
                   "gestational_age",
                   "birthweight",
                   "bmi_mum")
  #childhood
covariates_m10 <- c("smoking_mum",
                    "male_sex",
                    "education_mum",
                    "nulliparous",
                    "age_child",
                    "birthweight",
                    "bmi_mum")

#create SV function
SV.generate <- function(meth.matrix, pheno.data, variable.of.interest, model.covariates,n.sv, model_n){ 
  intersecting.samples <- intersect(pheno.data$Sample_Name, colnames(meth.matrix))
  pheno.data <- na.omit(pheno.data[which(pheno.data$Sample_Name %in% intersecting.samples), unique(c("Sample_Name", variable.of.interest, model.covariates))])
  meth.matrix <- meth.matrix[, match(pheno.data$Sample_Name, colnames(meth.matrix))]
  k = which(is.na(meth.matrix), arr.ind = TRUE)
  meth.matrix[k] = rowMedians(meth.matrix, na.rm = TRUE)[k[, 1]]
  mod = model.matrix(reformulate(paste0("pheno.data$",colnames(pheno.data[-1])))) #This eliminates samples.id
  #mod - FULL model matrix. The FULL model includes all the variables in the null model, as well as the variable of interest
  mod0 = mod[, -grep(paste0(variable.of.interest, collapse = "|"), colnames(mod))] #mod0 - NULL model matrix. The NULL model consist of known variables and covariates that must be included as adjustment variables
  sva.ret = sva(meth.matrix, mod = mod, mod0 = mod0, n.sv = n.sv)
  SVs = as.data.frame(sva.ret$sv)
  colnames(SVs) <-paste0("sv", 1:ncol(SVs), "_", model_n)
  cbind(Sample_Name = pheno.data$Sample_Name, SVs)
}

#add SVs to the phenotype dataset
pheno.sv.m10_Mono_child <- SV.generate(norm.beta.random, samplepheno[samplepheno$time_point==timepoint2,], "Mono_child", covariates_m10, n.sv=20, "m10_Mono_child")
pheno.sv.m10_CD8T_child <- SV.generate(norm.beta.random, samplepheno[samplepheno$time_point==timepoint2,], "CD8T_child", covariates_m10, n.sv=20, "m10_CD8T_child")
pheno.sv.m10_CD4T_child <- SV.generate(norm.beta.random, samplepheno[samplepheno$time_point==timepoint2,], "CD4T_child", covariates_m10, n.sv=20, "m10_CD4T_child")
pheno.sv.m10_NK_child <- SV.generate(norm.beta.random, samplepheno[samplepheno$time_point==timepoint2,], "NK_child", covariates_m10, n.sv=20, "m10_NK_child")
pheno.sv.m10_Bcell_child <- SV.generate(norm.beta.random, samplepheno[samplepheno$time_point==timepoint2,], "Bcell_child", covariates_m10, n.sv=20, "m10_Bcell_child")
pheno.sv.m10_Gran_child <- SV.generate(norm.beta.random, samplepheno[samplepheno$time_point==timepoint2,], "Gran_child", covariates_m10, n.sv=20, "m10_Gran_child")

sv_per_model <- mget(ls(pattern = "^pheno.sv.m"))
combined.svs <- reduce(sv_per_model, full_join, by = 'Sample_Name') 
names(combined.svs)
combined.pheno <- merge(samplepheno, combined.svs, by.x = "Sample_Name", by.y = "Sample_Name", all.x= T )
names(combined.pheno)

#---------------------------------------------------------------------#
#     MATCH DNA METHYLATION DATA TO PHENOTYPE DATA INCLUDING SVs      ----
#---------------------------------------------------------------------#

match_beta_func <- function(pheno_model) {
  #match betas to samplesheet
  m <- match(pheno_model$Sample_Name, colnames(norm.beta.random))
  #keep samples that match (this drops samples missing in phenotype that are available in betas (i.e., samples dropped when including SVs))
  norm.beta.random_model <- norm.beta.random[, m]
}

norm.beta.random_m10_Mono_child <- match_beta_func(pheno.sv.m10_Mono_child)
norm.beta.random_m10_CD8T_child <- match_beta_func(pheno.sv.m10_CD8T_child)
norm.beta.random_m10_CD4T_child <- match_beta_func(pheno.sv.m10_CD4T_child)
norm.beta.random_m10_NK_child <- match_beta_func(pheno.sv.m10_NK_child)
norm.beta.random_m10_Bcell_child <- match_beta_func(pheno.sv.m10_Bcell_child)
norm.beta.random_m10_Gran_child <- match_beta_func(pheno.sv.m10_Gran_child)

combined.beta <- match_beta_func(combined.pheno)

#check that rows in phenotype match columns in beta matrix
identical(pheno.sv.m10_Mono_child$Sample_Name, colnames(norm.beta.random_m10_Mono_child))
identical(pheno.sv.m10_CD8T_child$Sample_Name, colnames(norm.beta.random_m10_CD8T_child))
identical(pheno.sv.m10_CD4T_child$Sample_Name, colnames(norm.beta.random_m10_CD4T_child))
identical(pheno.sv.m10_NK_child$Sample_Name, colnames(norm.beta.random_m10_NK_child))
identical(pheno.sv.m10_Bcell_child$Sample_Name, colnames(norm.beta.random_m10_Bcell_child))
identical(pheno.sv.m10_Gran_child$Sample_Name, colnames(norm.beta.random_m10_Gran_child))


identical(combined.pheno$Sample_Name, colnames(combined.beta))


#check number of samples remaining in betas
dim(norm.beta.random_m10_Mono_child) 
dim(norm.beta.random_m10_CD8T_child) 
dim(norm.beta.random_m10_CD4T_child) 
dim(norm.beta.random_m10_NK_child) 
dim(norm.beta.random_m10_Bcell_child) 
dim(norm.beta.random_m10_Gran_child) 

dim(combined.beta) 

#---------------------------------------------------------------------#
#                        MODIFY VARIABLES FOR EWAS                    ----
#---------------------------------------------------------------------#

#change names of samplepheno and beta matrix
Pheno <- combined.pheno
meth <- combined.beta

# change variable names from the phenodata
Pheno$ID <- as.character(Pheno$Sample_Name)


#select cpg sites of interest
meth <- meth[c("cg19423175", "cg01500140", "cg13917614", "cg18247172", "cg22348356", "cg20674490"),]

# set the models

mod_form10a_cpg  <- formula("Mono_child ~ CpG + sv1_m10_Mono_child + sv2_m10_Mono_child + sv3_m10_Mono_child + sv4_m10_Mono_child + sv5_m10_Mono_child + sv6_m10_Mono_child + sv7_m10_Mono_child + sv8_m10_Mono_child + sv9_m10_Mono_child + sv10_m10_Mono_child + sv11_m10_Mono_child + sv12_m10_Mono_child + sv13_m10_Mono_child + sv14_m10_Mono_child + sv15_m10_Mono_child + sv16_m10_Mono_child + sv17_m10_Mono_child + sv18_m10_Mono_child + sv19_m10_Mono_child + sv20_m10_Mono_child + smoking_mum + male_sex + education_mum + nulliparous + age_child + birthweight + bmi_mum")
mod_form10b_cpg  <- formula("CD8T_child ~ CpG + sv1_m10_CD8T_child + sv2_m10_CD8T_child + sv3_m10_CD8T_child + sv4_m10_CD8T_child + sv5_m10_CD8T_child + sv6_m10_CD8T_child + sv7_m10_CD8T_child + sv8_m10_CD8T_child + sv9_m10_CD8T_child + sv10_m10_CD8T_child + sv11_m10_CD8T_child + sv12_m10_CD8T_child + sv13_m10_CD8T_child + sv14_m10_CD8T_child + sv15_m10_CD8T_child + sv16_m10_CD8T_child + sv17_m10_CD8T_child + sv18_m10_CD8T_child + sv19_m10_CD8T_child + sv20_m10_CD8T_child + smoking_mum + male_sex + education_mum + nulliparous + age_child + birthweight + bmi_mum")
mod_form10c_cpg  <- formula("CD4T_child ~ CpG + sv1_m10_CD4T_child + sv2_m10_CD4T_child + sv3_m10_CD4T_child + sv4_m10_CD4T_child + sv5_m10_CD4T_child + sv6_m10_CD4T_child + sv7_m10_CD4T_child + sv8_m10_CD4T_child + sv9_m10_CD4T_child + sv10_m10_CD4T_child + sv11_m10_CD4T_child + sv12_m10_CD4T_child + sv13_m10_CD4T_child + sv14_m10_CD4T_child + sv15_m10_CD4T_child + sv16_m10_CD4T_child + sv17_m10_CD4T_child + sv18_m10_CD4T_child + sv19_m10_CD4T_child + sv20_m10_CD4T_child + smoking_mum + male_sex + education_mum + nulliparous + age_child + birthweight + bmi_mum")
mod_form10d_cpg  <- formula("NK_child ~ CpG + sv1_m10_NK_child + sv2_m10_NK_child + sv3_m10_NK_child + sv4_m10_NK_child + sv5_m10_NK_child + sv6_m10_NK_child + sv7_m10_NK_child + sv8_m10_NK_child + sv9_m10_NK_child + sv10_m10_NK_child + sv11_m10_NK_child + sv12_m10_NK_child + sv13_m10_NK_child + sv14_m10_NK_child + sv15_m10_NK_child + sv16_m10_NK_child + sv17_m10_NK_child + sv18_m10_NK_child + sv19_m10_NK_child + sv20_m10_NK_child + smoking_mum + male_sex + education_mum + nulliparous + age_child + birthweight + bmi_mum")
mod_form10e_cpg  <- formula("Bcell_child ~ CpG + sv1_m10_Bcell_child + sv2_m10_Bcell_child + sv3_m10_Bcell_child + sv4_m10_Bcell_child + sv5_m10_Bcell_child + sv6_m10_Bcell_child + sv7_m10_Bcell_child + sv8_m10_Bcell_child + sv9_m10_Bcell_child + sv10_m10_Bcell_child + sv11_m10_Bcell_child + sv12_m10_Bcell_child + sv13_m10_Bcell_child + sv14_m10_Bcell_child + sv15_m10_Bcell_child + sv16_m10_Bcell_child + sv17_m10_Bcell_child + sv18_m10_Bcell_child + sv19_m10_Bcell_child + sv20_m10_Bcell_child + smoking_mum + male_sex + education_mum + nulliparous + age_child + birthweight + bmi_mum")
mod_form10f_cpg  <- formula("Gran_child ~ CpG + sv1_m10_Gran_child + sv2_m10_Gran_child + sv3_m10_Gran_child + sv4_m10_Gran_child + sv5_m10_Gran_child + sv6_m10_Gran_child + sv7_m10_Gran_child + sv8_m10_Gran_child + sv9_m10_Gran_child + sv10_m10_Gran_child + sv11_m10_Gran_child + sv12_m10_Gran_child + sv13_m10_Gran_child + sv14_m10_Gran_child + sv15_m10_Gran_child + sv16_m10_Gran_child + sv17_m10_Gran_child + sv18_m10_Gran_child + sv19_m10_Gran_child + sv20_m10_Gran_child + smoking_mum + male_sex + education_mum + nulliparous + age_child + birthweight + bmi_mum")


#################################################################
##########   NO NEED TO CHANGE ANYTHING AFTER THIS   ############
#################################################################


# make vector of all models
models <- ls()[grep("mod_form", ls())]

# save the models info
df <- data.frame(name = character(), formula = character())
for (i in 1:length(models)) {
  df[i, "name"] <- models[i]
  df[i, "formula"] <- deparse1(get(models[i]))
}

write.csv(df, file.path(workdir, paste0(cohort, "_models_", date, "_", analyst, ".csv")), row.names = FALSE)

# meth is a matrix of Illumina beta values for all samples(columns) and probes(rows)

# Match the methylation and phenotype data
meth <- meth[, na.omit(match(Pheno$ID, colnames(meth)))]
Pheno <- Pheno[match(colnames(meth), Pheno$ID),]
ifelse(all(Pheno$Row_names == colnames(meth)),
       "meth and phenotype data successfully matched",
       "Data not matched! Check what went wrong!")
Pheno <- droplevels(Pheno)

# Remove outliers using winsorization
winsorize <- function(methylation, pct = winsorize.pct) {
  quantiles <- matrixStats::rowQuantiles(methylation,
                                         probs = c(pct, 1-pct),
                                         na.rm = TRUE)
  low <- quantiles[,1]
  upper <- quantiles[,2]

  outliers.lower <- rowSums(methylation < low, na.rm = TRUE)
  outliers.upper <- rowSums(methylation > upper, na.rm = TRUE)

  idx <- which(methylation < low, arr.ind = TRUE)
  methylation[idx] <- low[idx[,1]]

  idx <- which(methylation > upper, arr.ind = TRUE)
  methylation[idx] <- upper[idx[,1]]

  n <- rowSums(!is.na(methylation))
  log <- data.frame(outliers.lower, outliers.upper, n)

  return(list(methylation = methylation, log = log))
}

# Replace outliers using winsorizing 1% (0.5% each side)
replace.outliers <- winsorize(meth, 0.005)
meth <- replace.outliers$methylation
outlier.log <- replace.outliers$log

# Save the log
write.csv(outlier.log,
          file = file.path(workdir, paste0(cohort, "_N_Winsorize0.01_Log_", date, "_", analyst, ".csv")),
          quote = FALSE,
          row.names = FALSE)

# Check to make sure the winsorization worked properly
str(meth)
nrow(meth) # number of lines (genes)
ncol(meth) # number of columns (conditions)
class(meth)
dim(meth)

# Check N cases with methylation data available (before excluding NA's)
# and order the files according to ID
index <- which(colnames(meth) %in% Pheno$ID)
length(index)
meth <- meth[, index]
meth <- meth[, order(colnames(meth))]
Pheno <- Pheno[order(Pheno$ID),]
ncol(meth)
nrow(meth)

# Export the distributions of methylation levels
meth_dist <- data.frame(Mean = rowMeans(meth, na.rm = TRUE),
                        St.Dev = apply(meth, 1, sd, na.rm = TRUE),
                        Min = apply(meth, 1, min, na.rm = TRUE),
                        Max = apply(meth, 1, max, na.rm = TRUE),
                        IQR = rowIQRs(meth, na.rm = TRUE),
                        row.names = rownames(meth))
filename <- paste0(cohort, "_methylation_distributions_", date, "_", analyst,".csv")
write.csv(meth_dist,
          file = file.path(workdir, filename),
          row.names = TRUE,
          quote = FALSE)

###                               ###
### Based on the PACE EWAS script ###
###                               ###

options(stringsAsFactors = FALSE, scipen = 999)
meth <- t(meth)
nprobes <- ncol(meth)

for(mod in models) {
  #Pre-allocate memory for the EWAS results
  results <- data.frame(probeid = colnames(meth),
                        beta = rep(0, times = nprobes),
                        se = rep(0, times = nprobes),
                        p_val = rep(0, times = nprobes),
                        z = rep(0, times = nprobes),
                        n = rep(0, times = nprobes))
  for(i in 1:nprobes) {
    tryCatch({
      CpG <- meth[, i]
      rlm.fit <- rlm(get(mod), data = Pheno, maxit = 200)
      test <- lmtest::coeftest(rlm.fit, sandwich::vcovHC(rlm.fit, type = "HC0"))
      
      BETA <- test[2, "Estimate"]
      SE <- test[2, "Std. Error"]
      PVAL <- test[2, "Pr(>|z|)"]
      N <- length(rlm.fit$residual)
      
      data.table::set(results, i, 2L, BETA)
      data.table::set(results, i, 3L, SE)
      data.table::set(results, i, 4L, PVAL)
      data.table::set(results, i, 5L, N)
    }, error = function(err) {
      message("Error in ", colnames(meth)[i])
      data.table::set(results, i, 2L, NA)
      data.table::set(results, i, 3L, NA)
      data.table::set(results, i, 4L, NA)
      data.table::set(results, i, 5L, NA)
    })
    
    if(i %% 4000 == 0) {cat("Progress:", 100*i/nprobes, "%\n")}
  }
  
  #Export results
  results <- results[order(results$p_val),]
  colnames(results) <- c("probeID", "BETA", "SE", "P_VAL", "N")
  mod_name <- sub("_form","", mod)
  filename <- paste0(cohort,"_",mod_name,"_",date,"_",analyst,".txt")
  write.table(results, 
              file = file.path(workdir, filename), 
              col.names = TRUE, 
              row.names = FALSE, 
              quote = FALSE)
  #gzip(filename)
  
  ##### Calculate lambda and sample size
  lambda <- qchisq(median(results$P_VAL, na.rm = T), df = 1, lower.tail = F) /
    qchisq(0.5, 1)
  N <- results$N[1]
  ## Save analysis sample sizes & lambda:
  samplesizelambda <- data.frame(cohort = character(), 
                                 model = character(), 
                                 lambda = numeric(), 
                                 N = numeric())
  samplesizelambda[1,] <- c(cohort, mod_name, lambda, N)
  write.table(samplesizelambda, 
              file.path(workdir, paste0(cohort, "_SampleSizeLambda_", mod_name,"_",date,"_",analyst, ".txt")), 
              quote = FALSE, 
              row.names = FALSE)
  

  
  
}
















