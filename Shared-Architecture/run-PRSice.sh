#!/bin/bash

#SBATCH --account=def-eparra
#SBATCH --time=10:00:00
#SBATCH --job-name=PRSice-tanning-PD
#SBATCH --output=PRSice-tanning-PD
#SBATCH --mail-user=c.abbatangelo@mail.utoronto.ca
#SBATCH --mail-type=ALL
#SBATCH --mem=20G

module load nixpkgs/16.09
module load gcc/7.3.0
module load r/3.6.1

R --file=PRSice_v1.25.R -q --args \
plink ./plink_1.9_linux_160914 \
base BASALCAR_DATA_PRS.txt \
target PD_DATA_PRS.txt \
size.targ 482730 \
covary F \
clump.snps F \
ld 1000kg_ref_eur \
best.thresh.on.bar F \
scatter.R2 T \
ggfig T \
barchart.levels 0.00000005,0.0000001,0.000001,0.00001,0.0001,0.001,0.05,0.1,0.3,0.5,1 \
bar.col.is.pval T \
figname BasalCellONPD \
slower 0.00000005 \
supper 1 \
sinc 0.005 \
fastscore T \
quantiles T \
quant.ref 1 \
num.quantiles 5 \
remove.mhc T \
report.best.score.only F \
report.individual.scores T \
no.regression T \
for.meta T \
sumsum T

#Format input files "as numeric"
#Note that in UKBB
  #ref: Reference allele on the forward strand.
  #alt: Alternate allele (not necessarily minor allele). Used as effect allele for GWAS.
  #So UKBB treats A2 as effect allele (Watanabe et al treats A1 as effect allele)
