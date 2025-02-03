### PERFORMING TWO-SAMPLE SUMMARY-DATA MR ANALYSIS WITH MRAPPS

##Step 0: Installation and loading packages

# 0.1. Add N to data if required:

Basal <- read_delim("BasalCellCarcinoma_RSID.txt", "\t", escape_double = FALSE, trim_ws = TRUE, progress = FALSE)
Hair <- read_delim("BrownHair_FullSumStats_WithBeta.txt", "\t", escape_double = FALSE, trim_ws = TRUE, progress = FALSE)
Melanoma <- read_delim("Melanoma_FullSumStats.txt", "\t", escape_double = FALSE, trim_ws = TRUE, progress = FALSE)
Nevi <- read_delim("Nevi_RSID.txt", "\t", escape_double = FALSE, trim_ws = TRUE, progress = FALSE)
RedHair <- read_delim("RedHair_FullSumStats_WithBeta.txt", "\t", escape_double = FALSE, trim_ws = TRUE, progress = FALSE)
Tanning <- read_delim("Tanning_FullSumStats.txt", "\t", escape_double = FALSE, trim_ws = TRUE, progress = FALSE)
VitD <- read_delim("VitaminDLevels_RSID.txt", "\t", escape_double = FALSE, trim_ws = TRUE, progress = FALSE)
Vitiligo <- read_delim("Vitiligo_FullSumStats_WithBeta.txt", "\t", escape_double = FALSE, trim_ws = TRUE, progress = FALSE)

head(Basal)
head(Hair)
head(Melanoma)
head(Nevi)
head(RedHair)
head(Tanning)
head(VitD)
head(Vitiligo)

Basal$N <- 456276
Hair$N <- 385603
Nevi$N <- 456348
RedHair$N <- 385603
Tanning$N <- 378364
VitD$N <- 417580
Vitiligo$N <- 44266

write_delim(Basal, "BasalCellCarcinoma_RSID_N.txt", delim = "\t", col_names = TRUE)
write_delim(Hair, "BrownHair_RSID_N.txt", delim = "\t", col_names = TRUE)
write_delim(Nevi, "Nevi_RSID_N.txt", delim = "\t", col_names = TRUE)
write_delim(RedHair, "RedHair_FullSumStats_N.txt", delim = "\t", col_names = TRUE)
write_delim(Tanning, "Tanning_FullSumStats_N.txt", delim = "\t", col_names = TRUE)
write_delim(VitD, "VitaminDLevels_RSID_N.txt", delim = "\t", col_names = TRUE)
write_delim(Vitiligo, "Vitiligo_FullSumStats_N.txt", delim = "\t", col_names = TRUE)

PD_data <- read_delim("PD_SumStats_RSID.txt", "\t", escape_double = FALSE, trim_ws = TRUE, progress = FALSE)

PD_data$N <- PD_data$N_controls + PD_data$N_cases

write_delim(PD_data, "PD_SumStats_RSID_N.txt", delim = "\t", col_names = TRUE)

# 0.2. Now proceed with installing MR-APPS:

install.packages("devtools")
install.packages("readr")
devtools::install_github("YangLabHKUST/MR-APSS")

library(devtools)
library(readr)
library(MRAPSS)

# You may have to remove and reinstall the package
#install.packages("cli")

# After installation, try loading the library and re-running your script
#library(cli)

## Step 1: Prepare data and estimate nuisance parameters

# 1.1. Format data

PD_raw <- readr::read_delim("PD_SumStats_RSID_N.txt", "\t",
                            escape_double = FALSE,
                            trim_ws = TRUE, progress = F)

Basal_raw <- readr::read_delim("BasalCellCarcinoma_RSID_N.txt", "\t", escape_double = FALSE,
                             trim_ws = TRUE, progress = F)

Hair_raw <- readr::read_delim("BrownHair_RSID_N.txt", "\t", escape_double = FALSE,
                               trim_ws = TRUE, progress = F)

Melanoma_raw <- readr::read_delim("Melanoma_FullSumStats.txt", "\t", escape_double = FALSE,
                              trim_ws = TRUE, progress = F)

Nevi_raw <- readr::read_delim("Nevi_RSID_N.txt", "\t", escape_double = FALSE,
                                  trim_ws = TRUE, progress = F)

RedHair_raw <- readr::read_delim("RedHair_FullSumStats_N.txt", "\t", escape_double = FALSE,
                              trim_ws = TRUE, progress = F)

Tanning_raw <- readr::read_delim("Tanning_FullSumStats_N.txt", "\t", escape_double = FALSE,
                                 trim_ws = TRUE, progress = F)

VitD_raw <- readr::read_delim("VitaminDLevels_RSID_N.txt", "\t", escape_double = FALSE,
                                 trim_ws = TRUE, progress = F)

Vitiligo_raw <- readr::read_delim("Vitiligo_FullSumStats_N.txt", "\t", escape_double = FALSE,
                              trim_ws = TRUE, progress = F)


head(PD_raw)
head(Basal_raw)
head(Hair_raw)
head(Melanoma_raw)
head(Nevi_raw)
head(RedHair_raw)
head(Tanning_raw)
head(VitD_raw)
head(Vitiligo_raw)

PD_form = format_data(PD_raw,
                      snp_col = "rsid",
                      b_col = "b",
                      se_col = "se",
                      freq_col = "freq",
                      A1_col = "A1",
                      A2_col = "A2",
                      p_col = "p",
                      n_col = "N")

Basal_form = format_data(Basal_raw,
                      snp_col = "RSID",
                      b_col = "B",
                      se_col = "SE",
                      freq_col = "FREQ",
                      A1_col = "A1",
                      A2_col = "A2",
                      p_col = "P",
                      n_col = "N")

Hair_form = format_data(Hair_raw,
                         snp_col = "SNPID_UKB",
                         b_col = "BETA",
                         se_col = "SE",
                         freq_col = "MAF",
                         A1_col = "A1",
                         A2_col = "A2",
                         p_col = "P",
                         n_col = "N")

Melanoma_form = format_data(Melanoma_raw,
                      snp_col = "variant_id",
                      b_col = "beta",
                      se_col = "standard_error",
                      freq_col = "effect_allele_frequency",
                      A1_col = "effect_allele",
                      A2_col = "other_allele",
                      p_col = "p_value",
                      n_col = "N")

Nevi_form = format_data(Nevi_raw,
                            snp_col = "RSID",
                            b_col = "B",
                            se_col = "SE",
                            freq_col = "FREQ",
                            A1_col = "A1",
                            A2_col = "A2",
                            p_col = "P",
                            n_col = "N")

RedHair_form = format_data(RedHair_raw,
                        snp_col = "SNPID_UKB",
                        b_col = "BETA",
                        se_col = "SE",
                        freq_col = "MAF",
                        A1_col = "A1",
                        A2_col = "A2",
                        p_col = "P",
                        n_col = "N")

Tanning_form = format_data(Tanning_raw,
                           snp_col = "SNPID_UKB",
                           b_col = "BETA",
                           se_col = "SE",
                           freq_col = "MAF",
                           A1_col = "A1",
                           A2_col = "A2",
                           p_col = "P",
                           n_col = "N")

VitD_form = format_data(VitD_raw,
                           snp_col = "RSID",
                           b_col = "B",
                           se_col = "SE",
                           freq_col = "FREQ",
                           A1_col = "A1",
                           A2_col = "A2",
                           p_col = "P",
                           n_col = "N")

Vitiligo_form = format_data(Vitiligo_raw,
                        snp_col = "SNP",
                        b_col = "BETA",
                        se_col = "SE",
                        freq_col = "MAF",
                        A1_col = "A1",
                        A2_col = "A2",
                        p_col = "P",
                        n_col = "N")

# 1.2. Harmonize the formatted data sets and estimate nuisance parameters

paras = est_paras(dat1 = PD_form,
                  dat2 = Nevi_form,
                  trait1.name = "PD",
                  trait2.name = "Nevi",
                  ldscore.dir = "./eur_w_ld_chr")

# 1.3. LD clumping

install.packages("ieugwasr")
library(ieugwasr)

options(ieugwasr_api = 'gwas-api.mrcieu.ac.uk/')

MRdat = clump(paras$dat,
              IV.Threshold = 0.15,
              SNP_col = "SNP",
              pval_col = "pval.exp",
              clump_kb = 1000,
              clump_r2 = 0.001)

## Step 2: Fit MR-APPS for causal inference

head(MRdat)
library(mvtnorm)

MRres = MRAPSS(MRdat,
               exposure = "Nevi",
               outcome = "PD",
               C = paras$C,
               Omega = paras$Omega,
               Cor.SelectionBias = T)

MRplot(MRres, exposure = "PD", outcome = "Vitiligo")

sensitivity(MRdat,
            Omega=paras$Omega,
            C = paras$C,
            exposure = "PD",
            outcome = "Vitiligo")
