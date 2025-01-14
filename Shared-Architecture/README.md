# Determination of Shared Genetic Architecture using Polygenic Risk Scores (PRS)

- A polygenic risk score (PRS) represents the cumulative impact of genetic variants across the genome, weighted by their effect sizes on a particular phenotype (Choi et al., 2020).
-  These effect sizes are typically derived from GWAS results, and only variants surpassing a specified p-value threshold are incorporated into the calculation.  The calculation of PRS often involves multiple thresholds (e.g., P-value thredhold = 5 x 10^-8, 1 x 10-5, 0.05, etc.) to ensure robust association detection. Originally applied to psychiatric conditions like schizophrenia and bipolar disorder (Eusden et al., 2015), this methodology has been adapted for other traits to disentangle complex genetic relationships.
-  Here we apply the first step of the program PRSice v1.25 (Eusden et al., 2015) to test for shared genetic architecture between PD and pigmentation traits.  This step does not calculate a polygenic risk score (PRS), rather it assesses how well a base trait can predict a target trait at different P-value thresholds.
-  Summary statistics for PD and the 9 pigmentation traits were used as input – the analysis was conducted bidirectionally, with PD serving as the base trait for each pigmentation phenotype and then as the target trait for each pigmentation phenotype.  SNPs in linkage disequielibrium were removed according to the PRSice default parameters and 10 principal components were used to control for population structure

## Additional resources

- Github: https://github.com/choishingwan/PRSice
- Website: https://choishingwan.github.io/PRSice/

## References

Choi, S. W., Mak, T. S. H., & O’Reilly, P. F. (2020). Tutorial: a guide to performing polygenic risk score analyses. Nature Protocols, 15(9). https://doi.org/10.1038/s41596-020-0353-1

Euesden, J., Lewis, C. M., & O'Reilly, P. F. (2016). PRSice: Polygenic Risk Score software v1.25. Bioinformatics.


