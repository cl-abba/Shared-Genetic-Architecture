install.packages("usethis")

# sessionInfo() make sure you are running on an R version 3+ and 64-bit platform
library(usethis)
usethis::edit_r_environ() #increase the physical an virtual memory so that R_MAX_VSIZE=51200 MB

# Install necessary pacakges
install.packages("remotes")
remotes::install_github("MRCIEU/TwoSampleMR")

# Library TwoSampleMR

library(TwoSampleMR)

redHair_input <- "RedHair_FullSumStats.txt"

# Read the .txt file into a data frame using read.table()
redHair_exposure_data <- read.table(redHair_input, header = TRUE, sep = "")

# Display the first few rows of the data frame
head(redHair_exposure_data)

    ## Add beta column if logistic 
    redHair_exposure_data$BETA = log(redHair_exposure_data$OR)
    head(redHair_exposure_data)
    
    write.table(file = 'RedHair_FullSumStats_WithBeta.txt', redHair_exposure_data, sep = "\t", quote = FALSE, row.names = FALSE)

# Filter based on pvalue threshold from PGS analysis 

redHair_exposure_data_pthresh <- subset(redHair_exposure_data, P <=0.002)

head(redHair_exposure_data_pthresh)

write.table(file = 'RedHair_Exposure_pthreshold0.002.txt', redHair_exposure_data_pthresh, sep = "\t", quote = FALSE, row.names = FALSE)

# Get instruments

exposure_redHair <- read_exposure_data(
  filename = "RedHair_Exposure_pthreshold0.002.txt",
  sep = "\t",
  snp_col = "SNPID_UKB",
  beta_col = "BETA",
  se_col = "SE",
  effect_allele_col = "A1",
  other_allele_col = "A2",
  eaf_col = "MAF_UKB",
  chr_col = "CHR",
  pval_col = "P",
  pos_col = "BP",
)

head(exposure_redHair)

# Clump

exposure_redHair_clumped <- clump_data(
  exposure_redHair,
  clump_kb = 10000,
  clump_r2 = 0.001,
  clump_p1 = 0.025,
  clump_p2 = 1,
  pop = "EUR"
)

write.table(exposure_redHair_clumped, file = "RedHair_Exposure_Clumped.txt", sep = "\t", row.names = FALSE, quote = FALSE)

PD_input <- "PD_FullSumStats_RSID.txt"

# Read the .txt file into a data frame using read.table()
PD_outcome_data <- read.table(PD_input, header = TRUE, sep = "")

head(PD_outcome_data)

# Get effects of instruments on outcome
outcome_PD <- read_outcome_data(
  snps = exposure_redHair_clumped$SNP,
  filename = "PD_FullSumStats_RSID.txt",
  sep = "\t",
  snp_col = "rsid",
  beta_col = "b",
  se_col = "se",
  effect_allele_col = "A1",
  other_allele_col = "A2",
  eaf_col = "freq",
  chr_col = "CHR",
  pval_col = "P",
  pos_col = "POS",
  ncase_col = "N_cases",
  ncontrol_col = "N_control"
)

head(outcome_PD)

# Harmonise the exposure and outcome data

redHair_on_PD_harmonized <- harmonise_data(
  exposure_dat = exposure_redHair_clumped, 
  outcome_dat = outcome_PD)

head(melanoma_on_PD_harmonized)

# Perform MR

results_redHair_on_PD <- mr(redHair_on_PD_harmonized)
write.table(file = 'results_redHair_on_PD.txt', results_redHair_on_PD, sep = "\t", quote = FALSE, row.names = FALSE)

het_redHair_on_PD <- mr_heterogeneity(redHair_on_PD_harmonized)
write.table(file = 'het_redHair_on_PD.txt', het_redHair_on_PD, sep = "\t", quote = FALSE, row.names = FALSE)

pleio_redHair_on_PD <- mr_pleiotropy_test(redHair_on_PD_harmonized)
write.table(file = 'pleio_redHair_on_PD.txt', pleio_redHair_on_PD, sep = "\t", quote = FALSE, row.names = FALSE)

ssnp_redHair_on_PD <- mr_singlesnp(redHair_on_PD_harmonized)
write.table(file = 'ssnp_redHair_on_PD.txt', ssnp_redHair_on_PD, sep = "\t", quote = FALSE, row.names = FALSE)

loo_redHair_on_PD <- mr_leaveoneout(redHair_on_PD_harmonized)
write.table(file = 'loo_redHair_on_PD.txt', loo_redHair_on_PD, sep = "\t", quote = FALSE, row.names = FALSE)

RAPS_redHair_on_PD <- mr(redHair_on_PD_harmonized, method_list = c("mr_raps"))

#### Remove outliers if there is evidence of pleiotropy or heterogeneity

# First ssnps file
head(ssnp_redHair_on_PD)

# Calculate the "z" column
ssnp_redHair_on_PD$z <- ssnp_redHair_on_PD$b / mean(ssnp_redHair_on_PD$b)

# Calculate the confidence interval for "z"
z_mean_ssnp <- mean(ssnp_redHair_on_PD$z)
z_sd_ssnp <- sd(ssnp_redHair_on_PD$z)
z_ci_upper_ssnp <- z_mean_ssnp + 1.96 * z_sd_ssnp
z_ci_lower_ssnp <- z_mean_ssnp - 1.96 * z_sd_ssnp

# Identify SNPs outside the 95% confidence interval
outliers_ssnp <- ssnp_redHair_on_PD$z < z_ci_lower_ssnp | ssnp_redHair_on_PD$z > z_ci_upper_ssnp
outlier_snps_ssnp <- ssnp_redHair_on_PD$SNP[outliers_ssnp]

# Print the list of SNPs outside the confidence interval
cat("SNPs outside the 95% confidence interval:\n")
cat(outlier_snps_ssnp, sep = "\n")

redHair_on_PD_harmonized_ssnps <- redHair_on_PD_harmonized[!redHair_on_PD_harmonized$SNP %in% outlier_snps_ssnp, ]

nrow(redHair_on_PD_harmonized_ssnps)

#### loo file now 

# Calculate the "z" column
loo_redHair_on_PD$z <- loo_redHair_on_PD$b / mean(loo_redHair_on_PD$b)

# Calculate the confidence interval for "z"
z_mean_loo <- mean(loo_redHair_on_PD$z)
z_sd_loo <- sd(loo_redHair_on_PD$z)
z_ci_upper_loo <- z_mean_loo + 1.96 * z_sd_loo
z_ci_lower_loo <- z_mean_loo - 1.96 * z_sd_loo

# Identify SNPs outside the 95% confidence interval
outliers_loo <- loo_redHair_on_PD$z < z_ci_lower_loo | loo_redHair_on_PD$z > z_ci_upper_loo
outlier_snps_loo <- loo_redHair_on_PD$SNP[outliers_loo]

# Print the list of SNPs outside the confidence interval
cat("SNPs outside the 95% confidence interval:\n")
cat(outlier_snps_loo, sep = "\n")

redHair_on_PD_harmonized_95CI <- redHair_on_PD_harmonized_ssnps[!redHair_on_PD_harmonized_ssnps$SNP %in% outlier_snps_loo, ]

nrow(redHair_on_PD_harmonized_95CI)

write.table(file = 'RedHairOnPD_Clumped_Harmonized_95CI.txt', redHair_on_PD_harmonized_95CI, sep = "\t", quote = FALSE, row.names = FALSE)

#### Re-test MR measurements

results_redHair_on_PD_95CI <- mr(redHair_on_PD_harmonized_95CI)
write.table(file = 'results_redHair_on_PD_95CI.txt', results_redHair_on_PD_95CI, sep = "\t", quote = FALSE, row.names = FALSE)

het_redHair_on_PD_95CI <- mr_heterogeneity(redHair_on_PD_harmonized_95CI)
write.table(file = 'het_redHair_on_PD_95CI.txt', het_redHair_on_PD_95CI, sep = "\t", quote = FALSE, row.names = FALSE)

pleio_redHair_on_PD_95CI <- mr_pleiotropy_test(redHair_on_PD_harmonized_95CI)
write.table(file = 'pleio_redHair_on_PD_95CI.txt', pleio_redHair_on_PD_95CI, sep = "\t", quote = FALSE, row.names = FALSE)

ssnp_redHair_on_PD_95CI <- mr_singlesnp(redHair_on_PD_harmonized_95CI)
write.table(file = 'ssnp_redHair_on_PD_95CI.txt', ssnp_redHair_on_PD_95CI, sep = "\t", quote = FALSE, row.names = FALSE)

loo_redHair_on_PD_95CI <- mr_leaveoneout(redHair_on_PD_harmonized_95CI)
write.table(file = 'loo_redHair_on_PD_95CI.txt', loo_redHair_on_PD_95CI, sep = "\t", quote = FALSE, row.names = FALSE)

RAPS_redHair_on_PD_95CI <- mr(redHair_on_PD_harmonized_95CI, method_list = c("mr_raps"))

#Plot the graphs (use CI if het was signfiant)

p1 <- mr_scatter_plot(results_melanoma_on_PD, melanoma_on_PD_harmonized)

p1[[1]]

length(p1)

p2 <- mr_forest_plot(single_snp_exp_red_hair_out_PD)
p2[[1]]
?mr_forest_plot

#Notes:
    #Clump in R first (BASE only)
    #Benefit of using MR for this is that there is an internal EUR LD referene, so we do not have to download from the cluster
    #Take list of SNPs from that output and use as input for --extract command in PRSice
    #Large window size and small r squared to ensure we are removing anything that is even remotely correlated with one another
    #No covariates
