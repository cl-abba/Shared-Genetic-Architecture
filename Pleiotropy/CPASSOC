# Set working directory
setwd("/Volumes/LaCie/MultiTrait")

# Load necessary libraries
library(ggplot2)
#install.packages("ggrepel")
library(ggrepel)
library(RColorBrewer)
library(dplyr)

# Define the list of corrected summary statistics files
input_files <- c(
  "basal_corrected_v2.txt",
  "brownhair_corrected_v2.txt",
  "melanoma_corrected_v2.txt",
  "nevi_corrected_v2.txt",
  "PD_corrected_v2.txt",
  "redhair_corrected_v2.txt",
  "tanning_corrected_v2.txt",
  "vitd_corrected_v2.txt",
  "vitiligo_corrected_v2.txt"
)

# Define the corresponding output file names for PRSice
output_files <- c(
  "basal_corrected_v2_prs.txt",
  "brownhair_corrected_v2_prs.txt",
  "melanoma_corrected_v2_prs.txt",
  "nevi_corrected_v2_prs.txt",
  "PD_corrected_v2_prs.txt",
  "redhair_corrected_v2_prs.txt",
  "tanning_corrected_v2_prs.txt",
  "vitd_corrected_v2_prs.txt",
  "vitiligo_corrected_v2_prs.txt"
)

# Loop through each input file, format for PRSice, and save the new files
for (i in seq_along(input_files)) {
  # Read the summary statistics file
  sum_stats <- read.table(input_files[i], header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  
  # Select and rename the relevant columns
  prs_data <- sum_stats[, c("SNP", "CHR", "BP", "A1", "A2", "B", "SE", "P")]
  colnames(prs_data) <- c("SNP", "CHR", "BP", "A1", "A2", "BETA", "SE", "P")
  
  # Write the reformatted data to the output file
  write.table(prs_data, output_files[i], quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
}

# Load custom functions for SHom and SHet
source("FunctionSet.R")

# Consolidate all summary statistics into one data frame 
files <- list.files(pattern = "_corrected_v2.txt") # List all processed summary statistics files
summary_data <- list()

for (i in seq_along(files)) {
  # Read each file
  df <- read.table(files[i], header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  
  # Rename Z column based on phenotype (Z1, Z2, ...)
  colnames(df)[colnames(df) == "Z"] <- paste0("Z", i)
  
  # Keep only relevant columns
  summary_data[[i]] <- df[, c("SNP", "CHR", "BP", paste0("Z", i))]
}

# Merge all data frames by SNP, retaining only SNPs common to all datasets
merged_data <- summary_data[[1]] # Start with the first dataset

for (i in 2:length(summary_data)) {
  # Perform an inner join to retain only SNPs present in all datasets
  merged_data <- merge(merged_data, summary_data[[i]][, c("SNP", paste0("Z", i))], by = "SNP", all = FALSE)
}

# Retain CHR and BP from the first dataset only
merged_data <- merged_data[, !duplicated(names(merged_data))]

# Write the consolidated data to a file
write.table(merged_data, "merged_summary_statistics_corrected_v2_Jan28.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

head(merged_data)
tail(merged_data)

# Load the pruned SNP list (no header)
pruned_snps <- read.table("pruned_EUR.prune.in", header = FALSE, stringsAsFactors = FALSE)

# Rename the column for clarity
colnames(pruned_snps) <- c("SNP")

# Filter merged_data to include only SNPs in the pruned list
filtered_data <- merged_data[merged_data$SNP %in% pruned_snps$SNP, ]

# Filter for Z-scores within [-1.96, 1.96]
filtered_data2 <- filtered_data[rowSums(sapply(filtered_data[, grep("^Z", names(filtered_data))], abs) <= 1.96) == ncol(filtered_data[, grep("^Z", names(filtered_data))]), ]

# Identify Z-score columns
z_columns <- grep("^Z", names(filtered_data2), value = TRUE)

# Create the correlation matrix
corMatrix <- diag(x = 1, ncol = length(z_columns), nrow = length(z_columns))

for (i in 1:(length(z_columns) - 1)) {
  for (j in (i + 1):length(z_columns)) {
    corMatrix[i, j] <- corMatrix[j, i] <- cor(filtered_data2[[z_columns[i]]], filtered_data2[[z_columns[j]]], use = "pairwise.complete.obs")
  }
}

# Prepare the full matrix for SHom and SHet tests
Sumstat <- merged_data[, z_columns]

# Sample sizes for traits
Samplesize <- c(456276, 385603, 456276, 456348, 482730, 385603, 378364, 417580, 44266)

# Perform SHom test
Test.shom <- SHom(Sumstat, Samplesize, corMatrix)

# P-values for SHom
p.shom <- pchisq(Test.shom, df = 1, ncp = 0, lower.tail = FALSE)
res.shom <- data.frame(SNP = merged_data$SNP, p.Shom = p.shom)

# Perform SHet test
# Estimate gamma distribution parameters
para <- EstimateGamma(N = 1E4, Samplesize, corMatrix)

# Split data into 20 folds for SHet testing
folds <- cut(seq(1, nrow(Sumstat)), breaks = 20, labels = FALSE)
res.shet <- list()

for (i in 1:20) {
  Indexes <- which(folds == i, arr.ind = TRUE)
  Test.shet <- SHet(Sumstat[Indexes, ], Samplesize, corMatrix)
  
  # P-values for SHet
  p.shet <- pgamma(q = Test.shet - para[3], shape = para[1], scale = para[2], lower.tail = FALSE)
  res.shet[[i]] <- data.frame(SNP = merged_data$SNP[Indexes], p.Shet = p.shet)
}

# Combine SHet results
res.shet <- do.call(rbind.data.frame, res.shet)

# Add chromosome and position to the results
head(res.shom)
res.shom <- merge(res.shom, merged_data[, c("SNP", "CHR", "BP")], by = "SNP")
res.shet <- merge(res.shet, merged_data[, c("SNP", "CHR", "BP")], by = "SNP")

# Save results
write.table(res.shom, "SHom_Results_corrected_v2_Jan28.txt", row.names = FALSE, sep = "\t", quote = FALSE)
write.table(res.shet, "SHet_Results_corrected_v2_Jan28.txt", row.names = FALSE, sep = "\t", quote = FALSE)

# Remember to change p.shet and p.shom to P for downstream analyses

library(qqman)       # For QQ and Manhattan plots
library(data.table)

data_het <- read.table("SHet_Results_corrected_v2_Jan28.txt", header = TRUE, stringsAsFactors = FALSE)
data_hom <- read.table("SHom_Results_corrected_v2_Jan28.txt", header = TRUE, stringsAsFactors = FALSE)

head(data_het)
head(data_hom)

manhattan(data_hom, main = "Manhattan Plot for SHom",
          genomewideline = -log10(5e-8), ylim = c(0, 80), cex = 0.45, cex.axis = 0.75,
          col = c("darkorchid4", "skyblue2"), chrlabs = c(1:22),
          xlab = "Chromosome", ylab = "-log10(p-value)")


# Filter data for p-values below the suggestive significance threshold (1e-5)
suggestive_hits_het <- data_het[data_het$P < 1e-5, ]

# Filter data for p-values below the genome-wide threshold (5e-8)
significant_hits_het <- data_het[data_het$P < 5e-8, ]

# Save the filtered data into a new text file
fwrite(suggestive_hits_het, "SHet_SuggestiveHits_AllPheno_corrected_v2_Jan28.txt", sep = "\t", quote = FALSE)
fwrite(significant_hits_het, "SHet_SignificantHits_AllPheno_corrected_v2_Jan28.txt", sep = "\t", quote = FALSE)

# Filter data for p-values below the suggestive significance threshold (1e-5)
suggestive_hits_hom <- data_hom[data_hom$P < 1e-5, ]

# Filter data for p-values below the genome-wide threshold (5e-8)
significant_hits_hom <- data_hom[data_hom$P < 5e-8, ]

# Save the filtered data into a new text file
fwrite(suggestive_hits_hom, "SHom_SuggestiveHitsAllPheno_corrected_v2_Jan28.txt", sep = "\t", quote = FALSE)
fwrite(significant_hits_hom, "SHom_SignficantHitsAllPheno_corrected_v2_Jan28.txt", sep = "\t", quote = FALSE)

### Generate merged file for comparing p values 

# Set working directory (adjust to your directory if needed)
#setwd("/path/to/your/directory")

# List of filenames corresponding to phenotypes
filenames <- c(
  "basal_corrected_v2.txt", 
  "brownhair_corrected_v2.txt", 
  "melanoma_corrected_v2.txt", 
  "nevi_corrected_v2.txt", 
  "pd_corrected_v2.txt", 
  "redhair_corrected_v2.txt", 
  "tanning_corrected_v2.txt", 
  "vitd_corrected_v2.txt", 
  "vitiligo_corrected_v2.txt"
)

# List of phenotype names corresponding to files
phenotypes <- c(
  "Basal", "BrownHair", "Melanoma", "Nevi", 
  "PD", "RedHair", "Tanning", "VitD", "Vitiligo"
)

# Initialize an empty list to store dataframes
data_list <- list()

# Load each file and extract required columns
for (i in seq_along(filenames)) {
  df <- read.table(filenames[i], header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  
  # Extract relevant columns and rename P and B columns with phenotype names
  df <- df[, c("SNP", "CHR", "BP", "A1", "A2", "P", "B")]
  colnames(df)[6:7] <- c(paste0(phenotypes[i], "_P"), paste0(phenotypes[i], "_B"))
  
  # Append to the list
  data_list[[i]] <- df
}

# Merge all dataframes by SNPID, allowing for missing values
#merged_data <- Reduce(function(x, y) merge(x, y, by = c("SNPID", "CHR", "POS"), all = TRUE), data_list)

merged_data <- Reduce(function(x, y) merge(x, y, by = c("SNP", "CHR", "BP", "A1", "A2"), all = TRUE), data_list)

# View the first and last rows of the merged dataframe
head(merged_data)
tail(merged_data)

# List of SNPs to filter (replace with relevant SNPs)
target_snps <- c(
  "rs79323519", "rs34789477", "rs10789406", "rs1172128", "rs4617998", 
  "rs10197246", "rs1432262", "rs13067032", "rs1155563", "rs356200", 
  "rs2331343", "rs6881722", "rs6895666", "rs9457478", "rs2737214", 
  "rs4644350", "rs10809826", "rs12419854", "rs4121403", "rs12821256", 
  "rs1279397", "rs9556380", "rs1885194", "rs72714109", "rs11632984", 
  "rs10459819", "rs72833470", "rs8077038", "rs2092378", "rs6059655"
)


# Filter the merged_data dataframe
filtered_data <- merged_data[merged_data$SNP %in% target_snps, ]

# View the filtered data
head(filtered_data)

# Save the filtered dataframe to a file
write.table(filtered_data, "SHom_TopSNPs_AllPheno_Jan28.txt", sep = "\t", row.names = FALSE, quote = FALSE)


# Save the merged dataframe to a .txt file
write.table(merged_data, "Merged_Phenotypes_corrected_v2_Feb3.txt", sep = "\t", row.names = FALSE, quote = FALSE)

### Identify pleiotropic loci according to CPASSOC only 

# Filter SHom significant SNPs
SHom_significant <- data_hom[data_hom$P < 5e-8, ]

# Filter SHet significant SNPs
SHet_significant <- data_hom[data_hom$P < 5e-8, ]

# Define regions for SHom
SHom_regions <- SHom_significant %>%
  mutate(start = BP - 500000, end = BP + 500000) %>%
  select(SNP, CHR, start, end)

# Define regions for SHet
SHet_regions <- SHet_significant %>%
  mutate(start = BP - 500000, end = BP + 500000) %>%
  select(SNP, CHR, start, end)

# Filter conventional GWAS SNPs with P < 5e-8 in any GWAS
Phenotypes_significant <- merged_data %>%
  filter_at(vars(ends_with("_P")), any_vars(. < 5e-8))

# Check overlaps for SHom
SHom_overlap <- Phenotypes_significant %>%
  inner_join(SHom_regions, by = c("CHR"), relationship = "many-to-many") %>%
  filter(BP >= start & BP <= end)

# Check overlaps for SHet
SHet_overlap <- Phenotypes_significant %>%
  inner_join(SHet_regions, by = c("CHR"), relationship = "many-to-many") %>%
  filter(BP >= start & BP <= end)

# Filter unique SHom SNPs
SHom_unique <- SHom_significant %>%
  filter(!(SNP %in% SHom_overlap$SNP))

# Filter unique SHet SNPs
SHet_unique <- SHet_significant %>%
  filter(!(SNP %in% SHet_overlap$SNP))

# Save SHom unique SNPs
write.table(SHom_unique, "SHom_Unique_Jan18.txt", row.names = FALSE, quote = FALSE, sep = "\t")

# Save SHet unique SNPs
write.table(SHet_unique, "SHet_Unique_Jan18.txt", row.names = FALSE, quote = FALSE, sep = "\t")

# Save SHom significant for Manhattan and downstream analyses
write.table(SHom_significant, "SHom_Significant_Jan20.txt", row.names = FALSE, quote = FALSE, sep = "\t")

# Save SHet significant for Manhattan labels and downstream analyses
write.table(SHet_significant, "SHet_Significant_Jan20.txt", row.names = FALSE, quote = FALSE, sep = "\t")

######Pairwise tests for pleiotropy

# Set working directory
setwd("/Volumes/LaCie/MultiTrait")

# Load necessary libraries
library(ggplot2)
#install.packages("ggrepel")
library(ggrepel)
library(RColorBrewer)
library(dplyr)

# Load custom functions for SHom and SHet
source("FunctionSet.R")

# Specify the files to be read
files <- c("basal_corrected_v2.txt", "PD_corrected_v2.txt") 

# Consolidate selected summary statistics into one data frame
summary_data <- list()

for (i in seq_along(files)) {
  # Read each specified file
  df <- read.table(files[i], header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  
  # Add Z-score column (B/SE)
  df$Z <- df$B / df$SE
  
  # Rename Z column based on phenotype (Z1, Z2, ...)
  colnames(df)[colnames(df) == "Z"] <- paste0("Z", i)
  
  # Keep only relevant columns
  summary_data[[i]] <- df[, c("SNP", "CHR", "BP", paste0("Z", i))]
}

# Merge the two data frames by SNPID, retaining only SNPs common to both datasets
merged_data <- merge(summary_data[[1]], summary_data[[2]][, c("SNP", "Z2")], by = "SNP", all = FALSE)

# Retain CHR and POS from the first dataset only
merged_data <- merged_data[, !duplicated(names(merged_data))]

head(merged_data)
tail(merged_data)

# Load the pruned SNP list (no header)
pruned_snps <- read.table("pruned_EUR.prune.in", header = FALSE, stringsAsFactors = FALSE)

# Rename the column for clarity
colnames(pruned_snps) <- c("SNP")

# Filter merged_data to include only SNPs in the pruned list
filtered_data <- merged_data[merged_data$SNP %in% pruned_snps$SNP, ]

# Filter for Z-scores within [-1.96, 1.96]
filtered_data2 <- filtered_data[rowSums(sapply(filtered_data[, grep("^Z", names(filtered_data))], abs) <= 1.96) == ncol(filtered_data[, grep("^Z", names(filtered_data))]), ]

# Identify Z-score columns
z_columns <- grep("^Z", names(filtered_data2), value = TRUE)

# Create the correlation matrix
corMatrix <- diag(x = 1, ncol = length(z_columns), nrow = length(z_columns))

for (i in 1:(length(z_columns) - 1)) {
  for (j in (i + 1):length(z_columns)) {
    corMatrix[i, j] <- corMatrix[j, i] <- cor(filtered_data2[[z_columns[i]]], filtered_data2[[z_columns[j]]], use = "pairwise.complete.obs")
  }
}

# Prepare the full matrix for SHom and SHet tests
Sumstat <- merged_data[, z_columns]

# Sample sizes for traits
Samplesize <- c(456276, 482730)

# Perform SHom test
Test.shom <- SHom(Sumstat, Samplesize, corMatrix)

# P-values for SHom
p.shom <- pchisq(Test.shom, df = 1, ncp = 0, lower.tail = FALSE)
res.shom <- data.frame(SNP = merged_data$SNP, p.Shom = p.shom)

# Perform SHet test
# Estimate gamma distribution parameters
para <- EstimateGamma(N = 1E4, Samplesize, corMatrix)

# Split data into 20 folds for SHet testing
folds <- cut(seq(1, nrow(Sumstat)), breaks = 20, labels = FALSE)
res.shet <- list()

for (i in 1:20) {
  Indexes <- which(folds == i, arr.ind = TRUE)
  Test.shet <- SHet(Sumstat[Indexes, ], Samplesize, corMatrix)
  
  # P-values for SHet
  p.shet <- pgamma(q = Test.shet - para[3], shape = para[1], scale = para[2], lower.tail = FALSE)
  res.shet[[i]] <- data.frame(SNP = merged_data$SNP[Indexes], p.Shet = p.shet)
}

# Combine SHet results
res.shet <- do.call(rbind.data.frame, res.shet)

# Add chromosome and position to the results
head(res.shom)
res.shom <- merge(res.shom, merged_data[, c("SNP", "CHR", "BP")], by = "SNP")
res.shet <- merge(res.shet, merged_data[, c("SNP", "CHR", "BP")], by = "SNP")

# Save results
write.table(res.shom, "SHom_Results_Basal_PD_Feb3.txt", row.names = FALSE, sep = "\t", quote = FALSE)
write.table(res.shet, "SHet_Results_Basal_PD_Feb3.txt", row.names = FALSE, sep = "\t", quote = FALSE)

# Remember to change p.shet and p.shom to P for downstream analyses

library(qqman)       # For QQ and Manhattan plots
library(data.table)

data_het <- read.table("SHet_Results_Basal_PD_Feb3.txt", header = TRUE, stringsAsFactors = FALSE)
data_hom <- read.table("SHom_Results_Basal_PD_Feb3.txt", header = TRUE, stringsAsFactors = FALSE)

head(data_het)
head(data_hom)

manhattan(data_het, main = "Pairwise Basal-PD Manhattan Plot for SHet",
          genomewideline = -log10(5e-8), ylim = c(0, 50), cex = 0.45, cex.axis = 0.75,
          col = c("seagreen3", "magenta4"), chrlabs = c(1:22),
          xlab = "Chromosome", ylab = "-log10(p-value)")

# Filter data for p-values below the suggestive significance threshold (1e-5)
suggestive_hits_het <- data_het[data_het$P < 1e-5, ]

# Filter data for p-values below the genome-wide threshold (5e-8)
significant_hits_het <- data_het[data_het$P < 5e-8, ]

# Save the filtered data into a new text file
fwrite(suggestive_hits_het, "SHet_SuggestiveHits_Basal_PD_Feb3.txt", sep = "\t", quote = FALSE)
fwrite(significant_hits_het, "SHet_SignificantHits_Basal_PD_Feb3.txt", sep = "\t", quote = FALSE)

# Filter data for p-values below the suggestive significance threshold (1e-5)
suggestive_hits_hom <- data_hom[data_hom$P < 1e-5, ]

# Filter data for p-values below the genome-wide threshold (5e-8)
significant_hits_hom <- data_hom[data_hom$P < 5e-8, ]

# Save the filtered data into a new text file
fwrite(suggestive_hits_hom, "SHom_SuggestiveHits_Basal_PD_Feb3.txt", sep = "\t", quote = FALSE)
fwrite(significant_hits_hom, "SHom_SignficantHits_Basal_PD_Feb3.txt", sep = "\t", quote = FALSE)

### Generate merged file for comparing p values 

# Set working directory (adjust to your directory if needed)
#setwd("/path/to/your/directory")

# List of filenames corresponding to phenotypes
filenames <- c(
  "basal_corrected_v2.txt", 
  "pd_corrected_v2.txt" 
)

# List of phenotype names corresponding to files
phenotypes <- c(
  "Basal", "PD"
)

# Initialize an empty list to store dataframes
data_list <- list()

# Load each file and extract required columns
for (i in seq_along(filenames)) {
  df <- read.table(filenames[i], header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  
  # Extract relevant columns and rename P and B columns with phenotype names
  df <- df[, c("SNP", "CHR", "BP", "A1", "A2", "P", "B")]
  colnames(df)[6:7] <- c(paste0(phenotypes[i], "_P"), paste0(phenotypes[i], "_B"))
  
  # Append to the list
  data_list[[i]] <- df
}

# Merge all dataframes by SNPI, allowing for missing values
merged_data <- Reduce(function(x, y) merge(x, y, by = c("SNP", "CHR", "BP", "A1", "A2"), all = TRUE), data_list)

# View the first and last rows of the merged dataframe
head(merged_data)
tail(merged_data)

# List of SNPs to filter (replace with relevant SNPs)
target_snps <- c(
  "rs145330152", "rs12075248", "rs4613239", "rs6719014", "rs11917652", "rs10513789", "rs34311866",
  "rs4389574", "rs356220", "rs9271377", "rs620490", "rs10905486", "rs112460025", "rs9897399", "rs214755"
)


# Filter the merged_data dataframe
filtered_data <- merged_data[merged_data$SNP %in% target_snps, ]

# View the filtered data
head(filtered_data)

# Save the filtered dataframe to a file
write.table(filtered_data, "SHom_TopSNPs_Basal_PD_Feb3.txt", sep = "\t", row.names = FALSE, quote = FALSE)

