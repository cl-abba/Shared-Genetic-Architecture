# Using Mendelian Randomization (MR) for causal inference analysis

## About MR 

- To assess possible causal effects between two traits of interest,  bidirectional MR analyses using STROBE guidelines can be applied.  For more on STROBE guidelines see Skrivankova et al., 2021 (https://pubmed.ncbi.nlm.nih.gov/34698778/).
- Two MR R packages are described in this repo:
    1) TwoSample MR v0.5.7 (https://mrcieu.github.io/TwoSampleMR/index.html) (Hemani et al., 2018)
    2) MR-APPS (https://github.com/YangLabHKUST/MR-APSS) (Hu et al., 2022).
- TwoSample MR computes MR results with the methods MR Egger, Weighted median, IVW, Simple mode and Weighted mode.
- MR-APPS is a novel MR method that accounts for pleiotropy, selection bias, population stratification and sample overlap.  In a recent benchmarking analysis, MR-APPS outperformed other MR methods, producing more accurate causal effect estimates with narrower confidence intervals (Hu et al., 2024).
- Exposure datasets in TwoSample MR were filtered based on the P-value threshold with the greatest variance explained (R2) from the corresponding PRS analysis, however, different methods of P-value thresholding can be implemented by the user.  Each exposure dataset was then clumped using 10,000kb windows with an r^2=0.001 to obtain independent instruments.  We used variant IDs, effect allele, other allele, effect allele frequency, effect size, standard error, and P-value for the analyses, as outlined by each program’s documentation.   Sensitivity analyses detecting the presence of heterogeneity and horizontal pleiotropy were performed using various test statistics in both packages.
- For more help with interpreting MR output, see the following helpful resource from Davies, Holmes & Davey Smith 2018 (https://pubmed.ncbi.nlm.nih.gov/30002074/)

## Additional support and tutorials

- Step by step tutorial for TwoSampleMR: https://mrcieu.github.io/TwoSampleMR/articles/introduction.html

## References

Hemani, G., Zheng, J., Elsworth, B., Wade, K. H., Haberland, V., Baird, D., et al. (2018). The MR-Base platform supports systematic causal inference across the human phenome. eLife, 7. https://doi.org/10.7554/eLife.34408
Hu, X., Cai, M., Xiao, J., Wan, X., Wang, Z., Zhao, H., & Yang, C. (2024). Benchmarking Mendelian Randomization methods for causal inference using genome-wide association study summary statistics. medRxiv. https://doi.org/10.1016/j.ajhg.2024.06.016
Hu, X., Zhao, J., Lin, Z., Wang, Y., Peng, H., Zhao, H., ... & Yang, C. (2022). Mendelian randomization for causal inference accounting for pleiotropy and sample structure using genome-wide summary statistics. Proceedings of the National Academy of Sciences, 119(28), e2106858119.




