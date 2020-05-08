# The proportion of variance associated by DNA methylation across traits in ALSPAC

## Broad steps:
1. Extract binary and continuous data from ALSPAC
2. Make methylation kinship matrices 
3. Use these to assess the proportion of variance associated with DNA methylation (h<sup>2</sup><sub>EWAS</sub>) across multiple traits using both the LDAK and GCTA models - assess difference between models
4. Examine the relationship between h<sup>2</sup><sub>EWAS</sub> and number of EWAS hits
5. Run sensitivity analyses to check some assumptions made when calculating h<sup>2</sup><sub>EWAS</sub> 

## Steps with code
1. Extract ALSPAC data [data_extraction.R](methyl_variance/R/data_extraction.R)
2. Check the distribution of the phenotypes and remove any phenotypes without normal distributions after rank normalisation [check_phenos.R](methyl_variance/R/check_phenos.R)
3. Filter CpG sites [filter_aries_cpgs.R](methyl_variance/R/filter_aries_cpgs.R)
4. Create files needed to produce kinship matrices [create_bim_fam_sp.R](methyl_variance/R/create_bim_fam_sp.R)
5. Generate top 10 genetic PCs [pca.sh](methyl_variance/Shell/pca.sh)
6. Gather covariates together and put them into the right format [Make_covar.R](methyl_variance/R/Make_covar.R)
7. Make methylation kinship matrices
	+ For GCTA model [make_kinship_TB.sh](methyl_variance/Shell/make_kinship_TB.sh), for LDAK model [make_kinship_mod2.sh](methyl_variance/Shell/make_kinship_mod2.sh), for checking influence of weighting CpGs by variance [make_kinship_sens.sh](methyl_variance/Shell/make_kinship_sens.sh)
8. Extract the results [extract_m2_estimates.R](methyl_variance/R/extract_m2_estimates.R)
9. Check kinship values - make new matrices with outliers removed for sensitivity analyses [Kinship_descriptives.R](methyl_variance/R/Kinship_descriptives.R)
10. Perform analysis to assess h<sup>2</sup><sub>EWAS</sub> across many traits in ALSPAC using both GCTA and LDAK models [GCTA_analysis.R](methyl_variance/R/GCTA_analysis.R)
11. Run sensitivity analyses for h<sup>2</sup><sub>EWAS</sub> estimates [GCTA_analysis.R](methyl_variance/R/GCTA_analysis.R)
12. Prune traits [prune_traits.R](methyl_variance/R/prune_traits.R)
13. Using [meffil](https://www.ncbi.nlm.nih.gov/pubmed/29931280), run EWAS on pruned traits [meffil_EWAS_script.r](EWAS/R/meffil_EWAS_script.r)
14. Extract DNAm sites that were identified using different P-value cutoff thresholds [extract_EWAS_hits.R](EWAS/R/extract_ewas_hits.R)
15. Assess differences between LDAK and GCTA models [gcta_ldak_comparison.R](methyl_variance/R/gcta_ldak_comparison.R)
16. Run the sensitivity analyses [sensitivity.R](methyl_variance/R/sensitivity.R)
17. Examine the relationship between h<sup>2</sup><sub>EWAS</sub> and number of EWAS hits - do this for both models! [post_EWAS.R](EWAS/R/post_EWAS.R)
	+ This involves just checking the association and assessing the ability of h<sup>2</sup><sub>EWAS</sub> to predict number of EWAS hits
