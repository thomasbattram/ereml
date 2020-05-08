# -----------------------------------------------------------------
# Checking the normality of phenotypes after rank transformation
# -----------------------------------------------------------------

pkgs <- c("tidyverse", "gridExtra", "readxl", "haven", "GenABEL")
lapply(pkgs, require, character.only=T)
devtools::load_all("~/repos/usefunc/")

FOM1_dat <- read_delim("methyl_variance/phen/ALSPAC/data_FOM.txt",
					   delim = "\t", locale = locale(encoding = "latin1"))
meta_dat <- read_tsv("methyl_variance/phen/ALSPAC/phenotype_metadata_FOM.txt")

# function to extract variable from the data, transform it and plot a histogram if continuous
plot_hist <- function(var) {
	print(var)
	x <- FOM1_dat[, var, drop = T]
	if (is.binary(x) | is.monomorphic(x)) {
		return(NULL)
	} else {
		print("transforming")
		x <- rntransform(x)
	}
	dat <- data.frame(x = x)
	plot <- ggplot(dat, aes(x = x)) +
		geom_histogram(bins = 30) +
		labs(title = var)
	return(plot)
}

not_int <- colnames(FOM1_dat)[map_lgl(FOM1_dat, is.character)]

## too many vars:
## checking whether shapiro wilk test can detect non-normality here or if it's likely too sensitive
shap_test <- lapply(FOM1_dat[,-grep(paste(not_int, collapse = "|"), colnames(FOM1_dat))], function(x) {
	if (!is.binary(x) & !is.monomorphic(x)) {
		shapiro.test(rntransform(x))
	}
})
shap_test2 <- shap_test[!map_lgl(shap_test, is.null)]

non_norm <- names(shap_test2)[map_lgl(shap_test2, function(x) x$p < 0.05/length(shap_test2))]

# Checking the distribution of the shapiro-wilk test p values and distribution of the traits
# as the shapiro-wilk test can often be too sensitive to any deviations 
# from normality with this many samples
pdf("pvals.pdf")
hist(pvals)
dev.off()
non_norm_dat <- FOM1_dat[, non_norm]

hist_list <- lapply(colnames(non_norm_dat), plot_hist)

pdf("checking_histograms.pdf")
marrangeGrob(hist_list, nrow=1, ncol=1)
dev.off()

# seems the shapiro-wilk test actually does a decent job with these rank transformed traits 
# so just writing out the phenotypes with evidence of non-normality
# according to the shapiro-wilk test  
write.table(non_norm, file = "methyl_variance/phen/ALSPAC/non_normal_traits.txt", quote = F, row.names = F, col.names = F, sep = "\t")