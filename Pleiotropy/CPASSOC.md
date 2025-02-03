# About CPASSOC

- We used the CPASSOC R package (version 1.01) (Zhu et al., 2015) to investigate single loci partly responsible for shared genetic architecture between PD and the nine pigmentation traits used in this study.  The CPASSOC package can be accessed at http://hal.case.edu/zhu-web/.
- To execute a CPASSOC analysis, a correlation matrix is necessary to adjust for phenotype correlations or those stemming from overlapping or related samples across different cohorts.  This matrix is estimated using summary statistics derived from independent SNPs in a genome-wide association study (GWAS).
- Zhu et al. (2015) recommend estimating these correlations using SNPs in linkage equilibrium. For datasets obtained from GWAS Catalog and GWAS Atlas, where only summary statistics are available, linkage disequilibrium (LD) patterns can be borrowed from external sources such as the 1000 Genomes Project (1KGP), which is available on the PLINK2 Resources page (https://www.cog-genomics.org/plink/2.0/resources#phase3_1kg).
- Following the methods outlined in Li & Zhu’s 2017 demonstration of CPASSOC, the SNP selection involved LD pruning at r2=0.2 using PLINK2 (Chang et al., 2015; Purcell & Chang, n.d.) (https://www.cog-genomics.org/plink/2.0/).
- SNPs with significant effects may bias correlations among summary statistics (so those with Z scores exceeding ±1.96) are excluded, following the methodology outlined by Li & Zhu (2017).
- CPASSOC calculates two different measures: SHom and SHet.  The SHom method is analogous to the fixed-effect meta-analysis approach (Willer, Li & Abecassis, 2010), but it incorporates adjustments for correlations in summary statistics across traits and cohorts, which may arise from related traits, overlapping datasets, or shared samples.
- SHet is an extension of SHom, permitting heterogeneity across trait effects.Given that CPASSOC combines multiple traits, the genome-wide significance threshold of P = 5 x 10-8 can be applied in SHom and SHet analyses as in a conventional GWAS to determine pleiotropic loci.
- The R code contains lines corresponding to the following:
-     Reading in the summary statistics
-     Creating the correlation matrix
-     Running the SHet and SHom tests (for all traits)
-     Creating Manhattan plots for SHet and SHom
-     Creating results files for SHet, SHom and lead SNPs
-     Repeating the procedure for pairwise tests for pleiotropy

## References

Chang, C. C., Chow, C. C., Tellier, L. C., Vattikuti, S., Purcell, S. M., & Lee, J. J. (2015). Second-generation PLINK: rising to the challenge of larger and richer datasets. Gigascience, 4(1), s13742-015.  https://doi.org/10.1186/s13742-015-0047-8

Li, X., & Zhu, X. (2017). Cross-phenotype association analysis using summary statistics from GWAS. Statistical Human Genetics: Methods and Protocols, 455-467.https://doi.org/ 10.1007/978-1-4939-7274-6_22

Purcell, S., & Chang, C. (n.d.). PLINK 2.0 [Computer software]. Retrieved January 14, 2025, from https://www.cog-genomics.org/plink/2.0/

Willer, C. J., Li, Y., & Abecasis, G. R. (2010). METAL: fast and efficient meta-analysis of genomewide association scans. Bioinformatics, 26(17), 2190-2191.  https://doi.org/10.1093/bioinformatics/btq340

Zhu, X., Feng, T., Tayo, B. O., Liang, J., Young, J. H., Franceschini, N., ... & Redline, S. (2015). Meta-analysis of correlated traits via summary statistics from GWASs with an application in hypertension. The American Journal of Human Genetics, 96(1), 21-36.  https://doi.ord/10.1016/j.ajhg.2014.11.011



