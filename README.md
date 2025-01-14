# Shared Genetic Architecture 
Scripts for investigating shared genetic architecuture (correlative relationships, causal inference and pleiotropy).

This repository describes steps in Mendelian Randomization (MR) analyses. Compliment to the paper, "Investigating shared genetic architecture between pigmentation genetics and Parkinson’s Disease" which represents the fifth chapter of my PhD thesis.

## Introduction

## Materials and methods

To replicate this protocol, the following materials are required:

1. Summary statistics for traits of interest
2. R/R studio, as well as storage and workspace that can be used to storing output files, running jobs and downloading programs

The methodology is as follows:

1. Genetic Correlation with LDSC (GenomicSEM)
2. Shared Architecture with Polygenic Risk Scores (PRS) (PRSice version 1.25)
3. Causal Inference with Mendelian Randomization (TwoSample MR and MR-APPS)
4. Pleiotropy/Multi-trait association analysis with CPASSOC (http://hal.case.edu/zhu-web/)
5. Gene set enrichment (FUMA https://fuma.ctglab.nl/ and ShinyGo http://bioinformatics.sdstate.edu/go74/)
