# -------------------------------------------------------
# Pruning traits
# -------------------------------------------------------

# Phenotypes in the analysis are correlated and this could bias 
# mean estimates of h2EWAS. This script assesses correlation between
# traits and extracts traits that have |r| < 0.2.

# -------------------------------------------------------
# setup
# -------------------------------------------------------

pkgs <- c("tidyverse", "GenABEL", "conflicted")
lapply(pkgs, require, character.only = T)
devtools::load_all("~/repos/usefunc")

setwd("methyl_variance")
dat <- read_tsv("phen/ALSPAC/data_FOM.txt")

# read in filtered likelihood estimates
lik_res <- read_tsv("results/FOM_filtered_m2_likelihood_estimates.txt")

blanket_phen <- unique(lik_res[lik_res$model == "blanket", "phen", drop = T])
grouping_phen <- unique(lik_res[lik_res$model == "grouping", "phen", drop = T])
overlap_phen <- intersect(blanket_phen, grouping_phen)

# Remove non normal traits as they're not going to be used in the analysis
non_norm <- read.delim("phen/ALSPAC/non_normal_traits.txt", header = F, stringsAsFactors=F)
dat <- dat %>%
	dplyr::select(-one_of(non_norm[[1]])) %>%
	dplyr::select(one_of(overlap_phen))

dim(dat)
# can remove or just keep them in for now and remove later if needed...
grep("questionnaire", colnames(dat), ignore.case = T, value = T)
# -------------------------------------------------------
# correlations
# -------------------------------------------------------

# remove monomorphic traits
mon_dat <- sapply(dat, is.monomorphic)
which(mon_dat)
dat <- dat[, !mon_dat]
dim(dat)
# seperate out binary and continuous phenotypes - easier to rank transform
is_bin <- sapply(dat, is.binary)
sum(is_bin)
bin_dat <- dat[, is_bin]
cont_dat <- dat[, !is_bin]

ids <- grep("aln|qlet", colnames(cont_dat), ignore.case = T, value = T)

# rank transform data (remove ids also)
rntdat <- apply(dplyr::select(cont_dat, -one_of(ids)), 2, rntransform)

comb_dat <- cbind(bin_dat, rntdat)

# Function to count the number of NAs in a matrix or data.frame
count_na <- function(dat, col_or_row = 2) {
	stopifnot(col_or_row %in% c(1,2))
	x <- apply(dat, col_or_row, function(x) {sum(is.na(x))})
	return(x)
}

# change things to factors
comb_dat[] <- lapply(seq_along(comb_dat), function(x) {
	out <- type.convert(comb_dat[[x]])
	out <- as.numeric(out)
	return(out)
})

cor_dat <- abs(cor(comb_dat, use = "pairwise.complete.obs")) # STILL LOADS OF NAs!!!!

naaa <- count_na(cor_dat)
sum(naaa > 0) # number phenotypes with at least one NA value in the correlation matrix: 2663
sum(naaa > 0.05 * nrow(cor_dat)) # removing those with greater than 5% NAs (arbitrary number)
rm_vars <- naaa > 0.05 * nrow(cor_dat)

cor_dat <- cor_dat[!rm_vars, !rm_vars]
dim(cor_dat)
naaa <- count_na(cor_dat)
sum(naaa > 0) # 2

# visualising the correlations! 
png("phenotype_correlation.png")
gplots::heatmap.2(cor_dat, trace = "none", scale = "none", cexRow = 0.001, dendrogram = "none")
dev.off()


# Extract traits that correlate with many others and remove them

# Function to extract the traits with a correlation above a certain value
extract_cor <- function(dat, cutoff = 0.4) {
	x <- apply(dat, 2, function(x) {sum(x > cutoff, na.rm = T)})
	y <- x[order(x, decreasing = T)]
	return(y)
}

cor_num <- extract_cor(cor_dat)

cor_num[1]
rm_list <- vector(mode = "character")
rm_num <- 1
temp_dat <- cor_dat
# loop takes the pheno correlated with the most variables removes it then starts again
while(cor_num[1] > 1) {
	print(cor_num[1])
	rm_list[[rm_num]] <- names(cor_num[1])
	rm_num <- rm_num + 1
	temp_dat <- temp_dat[!(rownames(temp_dat) %in% names(cor_num[1])), !(colnames(temp_dat) %in% names(cor_num[1]))]
	cor_num <- extract_cor(temp_dat)
}
length(colnames(temp_dat))

# here check if there seems to be any bias towards selecting phenotypes that are not
# correlated with another phenotype because of NAs...
naaa2 <- count_na(temp_dat)
summary(naaa2) # max is 12 NAs... and upper quartile is 0 NAs...
ind_var_cor_dat <- cor_dat[,colnames(temp_dat)]
naaa3 <- count_na(ind_var_cor_dat)
summary(naaa3) # upper quartile is 0...
summary(naaa) # not much difference between the two counts of NAs
# seems it doesn't make too much difference so just saving the results! 

save(temp_dat, file = "pruned_traits.RData")

# ph
png("phenotype_correlation_after_pruning.png")
gplots::heatmap.2(temp_dat, trace = "none", scale = "none", labRow = , dendrogram = "none")
dev.off()

### testing with ggplot!!! 
m <- temp_dat
new_dat <- data.frame(row=rownames(m)[row(m)], col=colnames(m)[col(m)], corr=c(m))

p <- ggplot(new_dat, aes(x = row, y = col, fill = corr)) +
	geom_tile() +
	theme(axis.text = element_blank())

ggsave("test_correlation_plot.png", plot = p)

m <- cor_dat
new_dat2 <- data.frame(row=rownames(m)[row(m)], col=colnames(m)[col(m)], corr=c(m))

p <- ggplot(new_dat2, aes(x = row, y = col, fill = corr)) +
	geom_tile() +
	theme(axis.text = element_blank())

ggsave("test_correlation_plot_all.png", plot = p)

