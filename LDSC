# Run LDSC using GenomicSEM in R

# LOAD PACKAGES ====
devtools::install_github("GenomicSEM/GenomicSEM")
library(GenomicSEM)
library(tidyverse)


# LOAD VARIANTS + SUMSTATS ====
variants <- read.delim("~/Desktop/WENDT_LAB/CRISTINA_MR/variants.tsv")

BasalCellCarcinoma_RSID <- read.delim("~/Desktop/WENDT_LAB/CRISTINA_MR/SUMSTATS/NEW_SUMSTATS/BasalCellCarcinoma_RSID.txt")
#BlondeHairSumStats <- read.delim("~/Desktop/WENDT_LAB/CRISTINA_MR/SUMSTATS/NEW_SUMSTATS/BlondeHairSumStats.txt")
BrownHair_FullSumStats_WithBeta <- read.delim("~/Desktop/WENDT_LAB/CRISTINA_MR/SUMSTATS/NEW_SUMSTATS/BrownHair_FullSumStats_WithBeta.txt")
#Eczema_FullSumStats <- read.delim("~/Desktop/WENDT_LAB/CRISTINA_MR/SUMSTATS/NEW_SUMSTATS/Eczema_FullSumStats.txt")
#EyeColour_FullSumStats <- read.delim("~/Desktop/WENDT_LAB/CRISTINA_MR/SUMSTATS/NEW_SUMSTATS/EyeColour_FullSumStats.txt")
Melanoma_FullSumStats <- read.delim("~/Desktop/WENDT_LAB/CRISTINA_MR/SUMSTATS/NEW_SUMSTATS/Melanoma_FullSumStats.txt")
Nevi_RSID <- read.delim("~/Desktop/WENDT_LAB/CRISTINA_MR/SUMSTATS/NEW_SUMSTATS/Nevi_RSID.txt")
PD_SumStats_RSID <- read.delim("~/Desktop/WENDT_LAB/CRISTINA_MR/SUMSTATS/NEW_SUMSTATS/PD_SumStats_RSID.txt")
RedHair_FullSumStats_WithBeta <- read.delim("~/Desktop/WENDT_LAB/CRISTINA_MR/SUMSTATS/NEW_SUMSTATS/RedHair_FullSumStats_WithBeta.txt")
SkinCol <- read.delim("~/Desktop/WENDT_LAB/CRISTINA_MR/SUMSTATS/NEW_SUMSTATS/SkinCol.txt")
Tanning_FullSumStats <- read.delim("~/Desktop/WENDT_LAB/CRISTINA_MR/SUMSTATS/NEW_SUMSTATS/Tanning_FullSumStats.txt")
VITD_DATA_PRS_2 <- read.delim("~/Desktop/WENDT_LAB/CRISTINA_MR/SUMSTATS/NEW_SUMSTATS/VITD_DATA_PRS_2.txt")
Vitiligo_FullSumStats_WithBeta <- read.delim("~/Desktop/WENDT_LAB/CRISTINA_MR/SUMSTATS/NEW_SUMSTATS/Vitiligo_FullSumStats_WithBeta.txt")


# FORMAT SUMSTATS + VARIANTS ====
# Use this to modify the SE when converting B to OR if necessary
#PD_subset <- PD_subset %>% mutate(OR = exp(b),
#                   SE = OR * se)

# VARIANTS
variants_subset <- variants %>% select(variant, chr, pos, ref, alt, rsid)
variants_eye_subset <- variants %>% select(rsid, alt)
#PD subset to merge with PD_subset
variants_pd_subset <- variants %>% select(rsid, chr, pos)
variants_pd_subset <- variants_pd_subset %>% unite(chr_pos, chr, pos)
colnames(variants_pd_subset)[1] <- "SNP"


# BASAL CELL CARCINOMA - LOGISTIC GWAS - n = 456,276
BasalCellCarcinoma_RSID$OR <- exp(BasalCellCarcinoma_RSID$B)
basal_subset <- BasalCellCarcinoma_RSID %>% select(RSID, A1, A2, OR, P, FREQ)
colnames(basal_subset)[6] <- "MAF"
basal_subset <- subset(basal_subset, basal_subset$RSID != is.na(basal_subset$RSID) &
                          basal_subset$A1 != is.na(basal_subset$A1) &
                          basal_subset$A2 != is.na(basal_subset$A2) &
                          basal_subset$OR != is.na(basal_subset$OR) & 
                          basal_subset$P != is.na(basal_subset$P) &
                          basal_subset$MAF != is.na(basal_subset$MAF))
write.table(basal_subset, file='basal_subset.txt', row.names = FALSE, col.names = TRUE, quote = F)


# BLONDE HAIR - LOGISTIC GWAS - n = 385,603
# NOT USING ANYMORE
# BlondeHairSumStats$OR <- exp(BlondeHairSumStats$Beta)
# blonde_subset <- BlondeHairSumStats %>% select(RSID, Effect_allele, Other_allele, OR, P_value)
# colnames(blonde_subset)[2] <- "A1"
# colnames(blonde_subset)[3] <- "A2"
# colnames(blonde_subset)[5] <- "P"
# blonde_subset <- subset(blonde_subset, blonde_subset$RSID != is.na(blonde_subset$RSID) &
#                           blonde_subset$A1 != is.na(blonde_subset$A1) &
#                           blonde_subset$A2 != is.na(blonde_subset$A2) &
#                           blonde_subset$OR != is.na(blonde_subset$OR) & 
#                           blonde_subset$P != is.na(blonde_subset$P))
# write.table(blonde_subset, file='blonde_subset.txt', row.names = FALSE, col.names = TRUE, quote = F)


# BROWN HAIR - LOGISTIC GWAS - n = 385,603
brown_subset <- BrownHair_FullSumStats_WithBeta %>% select(SNPID_UKB, A1, A2, OR, P, MAF)
colnames(brown_subset)[1] <- "rsID"
write.table(brown_subset, file='brown_subset.txt', row.names = FALSE, col.names = TRUE, quote = F)


# ECZEMA
#eczema_subset <- Eczema_FullSumStats %>% select(rsID, reference_allele, other_allele, beta, se, p.value)
#colnames(eczema_subset)[2] <- "A1"
#colnames(eczema_subset)[3] <- "A2"
#colnames(eczema_subset)[4] <- "B"
#colnames(eczema_subset)[5] <- "SE"
#colnames(eczema_subset)[6] <- "P"
# Following line removes any "INDEL" SNP IDs that are not following "rs#"
#eczema_subset <- eczema_subset[!grepl("INDEL", eczema_subset$rsID),]
#write.table(eczema_subset, file='eczema_subset.txt', row.names = FALSE, col.names = TRUE, quote = F)


# EYE COLOUR - LINEAR GWAS - n = 5,641
# NOT USING ANYMORE - SAMPLE SIZE TOO SMALL
# EyeColour_FullSumStats$SE <- EyeColour_FullSumStats$STD_FE/sqrt(5641)
# eyecolour_subset <- EyeColour_FullSumStats %>% select(RSID, A1, A2, BETA_FE, SE, PVALUE_FE)
# colnames(eyecolour_subset) <- c("rsID", "A1", "A2", "BETA", "SE", "P")
# write.table(eyecolour_subset, file='eyecolour_subset.txt', row.names = FALSE, col.names = TRUE, quote = F)


# MELANOMA - LOGISTIC GWAS - n = 456,348
Melanoma_FullSumStats$OR <- exp(Melanoma_FullSumStats$beta)
Melanoma_FullSumStats$MAF <- 1 - Melanoma_FullSumStats$effect_allele_frequency
melanoma_subset <- Melanoma_FullSumStats %>% select(variant_id, effect_allele, other_allele, OR, p_value, MAF, N)
colnames(melanoma_subset)[1] <- "rsID"
colnames(melanoma_subset)[2] <- "A1"
colnames(melanoma_subset)[3] <- "A2"
colnames(melanoma_subset)[5] <- "P"
write.table(melanoma_subset, file='melanoma_subset.txt', row.names = FALSE, col.names = TRUE, quote = F)


# NEVI - LOGISTIC GWAS - n = 456,348
Nevi_RSID$OR <- exp(Nevi_RSID$B)
Nevi_RSID$MAF <- 1 - Nevi_RSID$FREQ
nevi_subset <- Nevi_RSID %>% select(RSID, A1, A2, OR, P, FREQ)
colnames(nevi_subset)[6] <- "MAF"
write.table(nevi_subset, file='nevi_subset.txt', row.names = FALSE, col.names = TRUE, quote = F)

# PD - LOGISTIC GWAS - n = 482,720
PD_SumStats_RSID$OR <- exp(PD_SumStats_RSID$b)
PD_SumStats_RSID$N <- PD_SumStats_RSID$N_cases + PD_SumStats_RSID$N_controls
PD_SumStats_RSID$MAF <- 1 - PD_SumStats_RSID$freq
PD_subset <- PD_SumStats_RSID %>% select(rsid, A1, A2, OR, p, N, MAF)
#PD_subset$SNP<-gsub("chr","",as.character(PD_subset$SNP))
#separate(data = PD_subset, col = SNP, into = c("CHR", "POS"), sep = "_")
#PD_subset <- merge(PD_subset, variants_pd_subset, by = "SNP")
write.table(PD_subset, file='PD_subset.txt', row.names = FALSE, col.names = TRUE, quote = F)


# RED HAIR - LOGISTIC GWAS - n = 385,603
red_subset <- RedHair_FullSumStats_WithBeta %>% select(SNPID_UKB, A1, A2, OR, P, MAF)
colnames(red_subset)[1] <- "rsID"
write.table(red_subset, file='red_subset.txt', row.names = FALSE, col.names = TRUE, quote = F)


# SKIN COLOUR - LINEAR GWAS - n = 381,433
skin_subset <- SkinCol %>% select(SNPID_UKB, A1, A2, BETA, SE, P, MAF)
colnames(skin_subset)[1] <- "rsID"
write.table(skin_subset, file='skin_subset.txt', row.names = FALSE, col.names = TRUE, quote = F)


# TANNING - LINEAR GWAS - n = 378,364
tanning_subset <- Tanning_FullSumStats %>% select(SNPID_UKB, A1, A2, BETA, SE, P, MAF)
colnames(tanning_subset)[1] <- "rsID"

write.table(tanning_subset, file='tanning_subset.txt', row.names = FALSE, col.names = TRUE, quote = F)

# Single lambda correction for Tanning SE using methods from Georgiopoulos et al. 2016 (doi.org/10.1017/S0016672316000069)
# Run LDSC first and make sure the lambda value is correct as below, if >> 1
sqrt_lambda <- sqrt(2.2174)
tanning_subset$SE <- tanning_subset$SE*sqrt_lambda


# VITAMIN D - LINEAR GWAS - n = 417,580
VITD_DATA_PRS_2$MAF <- 1 - VITD_DATA_PRS_2$FREQ
vitd_subset <- VITD_DATA_PRS_2 %>% select(SNP, A1, A2, BETA, SE, P, MAF)
colnames(vitd_subset)[7] <- "MAF"
write.table(vitd_subset, file='vitd_subset.txt', row.names = FALSE, col.names = TRUE, quote = F)


# VITILIGO - LINEAR GWAS - n = 44,266
vitiligo_subset <- Vitiligo_FullSumStats_WithBeta %>% select(SNP, A1, A2, BETA, SE, P, MAF)
colnames(vitiligo_subset)[1] <- "rsID"
write.table(vitiligo_subset, file='vitiligo_subset.txt', row.names = FALSE, col.names = TRUE, quote = F)


## MUNGE STEP ====
# This munge step is identical to the one outlined in the gSEM github page

# LOAD NEW TXT FILES
basal_subset <- read.csv("~/Desktop/WENDT_LAB/CRISTINA_MR/BN_LDSC_INPUTS/basal_subset.txt", sep="")
#blonde_subset <- read.csv("~/Desktop/WENDT_LAB/CRISTINA_MR/BN_LDSC_INPUTS/blonde_subset.txt", sep="")
brown_subset <- read.csv("~/Desktop/WENDT_LAB/CRISTINA_MR/BN_LDSC_INPUTS/brown_subset.txt", sep="")
#eyecolour_subset <- read.csv("~/Desktop/WENDT_LAB/CRISTINA_MR/BN_LDSC_INPUTS/eyecolour_subset.txt", sep="")
melanoma_subset <- read.csv("~/Desktop/WENDT_LAB/CRISTINA_MR/BN_LDSC_INPUTS/melanoma_subset.txt", sep="")
nevi_subset <- read.csv("~/Desktop/WENDT_LAB/CRISTINA_MR/BN_LDSC_INPUTS/nevi_subset.txt", sep="")
PD_subset <- read.csv("~/Desktop/WENDT_LAB/CRISTINA_MR/BN_LDSC_INPUTS/PD_subset.txt", sep="")
red_subset <- read.csv("~/Desktop/WENDT_LAB/CRISTINA_MR/BN_LDSC_INPUTS/red_subset.txt", sep="")
skin_subset <- read.csv("~/Desktop/WENDT_LAB/CRISTINA_MR/BN_LDSC_INPUTS/skin_subset.txt", sep="")
tanning_subset <- read.csv("~/Desktop/WENDT_LAB/CRISTINA_MR/BN_LDSC_INPUTS/tanning_subset.txt", sep="")
vitd_subset <- read.csv("~/Desktop/WENDT_LAB/CRISTINA_MR/BN_LDSC_INPUTS/vitd_subset.txt", sep="")
vitiligo_subset <- read.csv("~/Desktop/WENDT_LAB/CRISTINA_MR/BN_LDSC_INPUTS/vitiligo_subset.txt", sep="")

# Create vector of the summary statistics files
files <- c("basal_subset.txt", "brown_subset.txt", 
           "melanoma_subset.txt", "nevi_subset.txt", 
           "PD_subset.txt", "red_subset.txt", "skin_subset.txt", 
           "tanning_subset.txt", "vitd_subset.txt", "vitiligo_subset.txt")

# Define the reference file being used to allign alleles across summary stats
#Here we are using hapmap3
hm3 <- "eur_w_ld_chr/w_hm3.snplist"

# Name the traits 
trait.names <- c("basal", "brown", 
                 "melanoma", "nevi",
                 "PD", "red", "skin", 
                 "tanning", "vitd", "vitiligo")

# List the sample sizes. If SNP-specific sum of effective sample sizes, use NA
# Melanoma and PD have SNP-specific sample sizes
N = c(456276, 385603, 
      NA, 456348, 
      NA, 385603, 381433, 
      378364, 417580, 44266)

# Define imputation quality filter
info.filter = 0.9

# Define MAF filter
maf.filter = 0.01

# Munge
munge(files = files, hm3 = hm3, trait.names = trait.names, N = N, info.filter = info.filter, maf.filter = maf.filter)


## LDSC =====
# This LDSC step is identical to the one outlined in the gSEM github page

# Vector of munged summary statisitcs
traits <- c("basal.sumstats.gz", "brown.sumstats.gz", 
            "melanoma.sumstats.gz", "nevi.sumstats.gz", 
            "PD.sumstats.gz", "red.sumstats.gz", "skin.sumstats.gz", 
            "tanning.sumstats.gz", "vitd.sumstats.gz", "vitiligo.sumstats.gz")

# Enter sample prevalence of .5 to reflect that all traits were munged using the sum of effective sample size
# If traits are continuous (linear GWAS model), use NA
sample.prev <- c(4257/456276, 146394/385603, 
                 3564/456348, 351/456348, 
                 33674/482730, 17329/385603, NA, 
                 NA, NA, 4680/44266)

#vector of population prevalences
population.prev <- c(0.00089,0.3103,
                     0.135, 0.015,
                     0.00572, 0.04495, NA, 
                     NA, NA, 0.008)

#the folder of LD scores
ld <- "eur_w_ld_chr/"

#the folder of LD weights [typically the same as folder of LD scores]
wld <- "eur_w_ld_chr/"

#name the traits
trait.names <- c("basal", "brown", 
                "melanoma", "nevi",
                 "PD", "red", "skin", 
                 "tanning", "vitd", "vitiligo")

#run LDSC
LDSCoutput <- ldsc(traits = traits, sample.prev = sample.prev, population.prev = population.prev, ld = ld, wld = wld, trait.names = trait.names)

#optional command to save the output as a .RData file for later use
save(LDSCoutput, file="20240307_BN_CRISTINAMR_LDSCOutput.RData")

#Lambda for tanning summary statistics was very high (greater than 2)
