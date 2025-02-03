#!/bin/bash

#SBATCH --account=def-eparra
#SBATCH --time=5:00:00
#SBATCH --job-name=prep_ref
#SBATCH --output=prep_ref.out
#SBACTH --error=prep_ref.err
#SBATCH --mail-user=c.abbatangelo@mail.utoronto.ca
#SBATCH --mail-type=ALL
#SBATCH --mem=20G

module load nixpkgs/16.09
module load plink/2.00-10252019-avx2

# Specify input VCF file
input_vcf="g1000_eur.vcf.gz"

# Specify output base file name (without extensions)
output_base="1000kg-ref-eur"

# Convert VCF to PLINK binary format
plink2 --vcf "$input_vcf" --make-bed --out "$output_base"


