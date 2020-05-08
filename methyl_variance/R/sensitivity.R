# -------------------------------------------------------
# Sensitivity analyses
# -------------------------------------------------------

#
# Purpose of script: Examine how changing the methylation kinship matrix by: 
# 1. Changing outlier cutoff and 2. changing CpG weighting alters m2 estimates
# and overall trend for m2 estimates to be associated with number of EWAS hits
#

pkgs <- c("tidyverse", "cowplot", "gridExtra")
lapply(pkgs, require, character.only = T)

devtools::load_all("~/repos/usefunc/")

setwd("methyl_variance")
timepoints="FOM"

# Read in phenotypes extracted from ALSPAC
alsp_phen <- read_tsv("phen/ALSPAC/data_FOM.txt")
phens <- colnames(alsp_phen)[-grep("aln|qlet", colnames(alsp_phen))]

# remove non-normal phenotypes
non_norm <- read.delim("phen/ALSPAC/non_normal_traits.txt", header = F, stringsAsFactors=F)
non_norm <- non_norm[[1]]
phens <- phens[!phens %in% non_norm]

# Independent vars only!!!
load("pruned_traits.RData")
ind_vars <- rownames(temp_dat)

# -------------------------------------------------------
# 1. Altering methylation kinship matrix diagonal outlier cutoff
# -------------------------------------------------------

# Extract ldak outputs
path <- paste0("ldak_estimates/", timepoints, "/sens/")
files <- list.files(path = path)
files <- grep("cutoff", files, value = T)
files <- grep("reml", files, value = T)
length(files)
dat_name <- paste0("results/sens/FOM_outlier_ldak_output.RData")

parameters <- tibble(file = files) %>%
	mutate(phen = gsub(".*cutoff_", "", files)) %>%
	mutate(phen = gsub("\\.reml", "", phen)) %>%
	dplyr::filter(phen %in% ind_vars) # only independent variables! 

# Extract phenotype and cutoff from file names
out_cut <- gsub("cutoff_.*", "", parameters$file)
out_cut <- gsub(".*_", "", out_cut)
unique_nam <- paste(parameters$phen, out_cut, sep = "_CUT_")

# check that the separater isn't in the phenotype names
grep("_CUT_", unique_nam, value = T)

out_res <- lapply(parameters$file, function(file) {	
	# Read in results file
	dat <- read_ldak(paste0(path, file))
	return(dat)
})


# name the results
names(out_res) <- unique_nam

# add the original results that did not remove any individuals! 
load("results/FOM_inf_mod_ldak_output.RData")
res <- res[names(res) %in% ind_vars]
names(res) <- paste0(names(res), "_CUT_NA")
head(res)
# Combine results
ldak_res <- c(res, out_res)

save(ldak_res, file = dat_name)

# load in the data
load(dat_name)

# extract the likelihood data
lik_res <- ldak_res %>%
	map_df("likelihood_estimates", .id = "unique_id") %>%
	separate(col = unique_id, into = c("phen", "cut"), sep = "_CUT_") %>%
	dplyr::filter(variable == "Alt_Likelihood") %>%
	mutate(cut = ifelse(cut == 0, NA, cut))

# summarise likelihoods in each case
lik_res %>%
	group_by(cut) %>%
	summarise(med_lik = median(value))

# Count number of results when using different outlier cutoffs
cut_count <- lik_res %>% count(cut)
phen_count <- lik_res %>% count(phen)

# Make some plots
p1_lik <- ggplot(lik_res, aes(x = cut)) +
	geom_histogram(stat = "count") +
	labs(title = "Successful REML analyses at differing methylation kinship matrix diagonal thresholds",
		x = "Outlier threshold cutoff")

ggsave("results/sens/outlier_cutoff_hist.pdf", plot = p1_lik)

# Extract the phenotypes where the analysis ran for all matrices
phen_out <- phen_count[phen_count$n == 8, "phen"][[1]]
res_out <- filter(lik_res, phen %in% phen_out)
dim(res_out)
summary(res_out$value)
# crazy likehood value --> removing
res_out <- res_out %>%
	dplyr::filter(value < max(value))

# likelihood varying with cut
p2_lik <- ggplot(res_out, aes(x = cut, y = value)) +
	geom_boxplot() +
	labs(x = "Outlier threshold cutoff", 
		 y = "likelihood") +
	theme_bw()

ggsave("results/sens/outlier_cutoff_likelihood_boxplot.pdf", plot = p2_lik)


# extract the m2 values so we can compare m2 with cut
m2_res <- ldak_res %>%
	map_df("her_estimates", .id = "unique_id") %>%
	separate(col = unique_id, into = c("phen", "cut"), sep = "_CUT_") %>%
	dplyr::filter(Component == "Her_ALL") %>%
	mutate(cut = ifelse(cut == 0, NA, cut)) %>%
	rename(m2 = Heritability, se = Her_SD) %>%
	dplyr::select(phen, m2, se, cut)



# m2 varying with cut
p1_m2 <- ggplot(m2_res, aes(x = cut, y = m2)) +
	geom_boxplot() +
	labs(title = "Association between differing methylation kinship matrix diagonal thresholds and m2",
		x = "Outlier threshold cutoff", 
		y = bquote(m^2)) +
	theme_bw()

ggsave("results/sens/outlier_cutoff_m2_boxplot.pdf", plot = p1_m2)

temp_res <- dplyr::filter(m2_res, m2 < 0.3) %>%
	dplyr::filter(m2 > -0.3)

p1.5_m2 <- ggplot(temp_res, aes(x = cut, y = m2)) +
	geom_boxplot() +
	labs(title = "Association between differing methylation kinship matrix diagonal thresholds and m2",
		x = "Outlier threshold cutoff", 
		y = bquote(m^2))

ggsave("results/sens/outlier_cutoff_m2_boxplot_zoomed.pdf", plot = p1.5_m2)

set.seed(2)
#phen_out <- phen_count[phen_count$nn == 8, "phen"]
res_out <- filter(m2_res, phen %in% sample(phen, 10))

p2_m2 <- ggplot(res_out, aes(x = cut, y = m2, group = phen, colour = phen)) +
	geom_line() +
	geom_point() +
	theme(legend.position = "none") +
	labs(title = "Association between differing methylation kinship matrix diagonal thresholds and m2",
		x = "Outlier threshold cutoff", 
		y = bquote(m^2))

ggsave("results/sens/outlier_cutoff_m2_line.pdf", plot = p2_m2)

# Look at association of cutoff threshold and trend (assoc of m2 and EWAS hit number) 

# Read in some EWAS hits data at e-5
ewas_dat <- read.delim(paste0("../EWAS/results/", timepoints, "/derived/top_hits_e-5.txt"), stringsAsFactors = F)

# extract hits
phen_hits <- as.data.frame(table(ewas_dat$phen)) %>%
	rename(phen = Var1, hit_count = Freq) %>%
	dplyr::filter(phen %in% ind_vars)

# join hits to ldak res
hit_res <- left_join(phen_hits, m2_res) %>%
	mutate(p = pnorm(-abs(m2 / se))) %>% 
	mutate(p = ifelse(m2 < 0, 1, p)) %>%
	mutate(var_exp = ifelse(p < 0.05, "Yes", "No"))
head(hit_res)

p3_m2 <- ggplot(hit_res, aes(x = var_exp, y = log10(hit_count))) +
	geom_boxplot() +
	facet_grid(. ~ cut) +
	labs(x = bquote(m^2 ~ " > 0 (P < 0.05)"))

ggsave("results/sens/outlier_cutoff_hit_count_box.pdf", plot = p3_m2)

# -------------------------------------------------------
# 2. Altering methylation kinship matrix alpha value (weighting of predictors)
# -------------------------------------------------------

# Extract ldak outputs
path <- paste0("ldak_estimates/", timepoints, "/sens/")
files <- list.files(path = path)
files <- grep("alpha", files, value = T)
files <- grep("reml", files, value = T)
head(files)
length(files)
dat_name <- paste0("results/sens/", timepoints, "_alpha_ldak_output.RData")

parameters <- tibble(file = files) %>%
	mutate(phen = gsub(".*alpha_", "", files)) %>%
	mutate(phen = gsub("\\.reml", "", phen)) %>%
	dplyr::filter(phen %in% ind_vars) # only independent variables! 

# Extract phenotype and cutoff from file names
out_alpha <- gsub("_alpha_.*", "", parameters$file)
out_alpha <- gsub(".*_", "", out_alpha)
unique_nam <- paste(parameters$phen, out_alpha, sep = "_ALPHA_")

# check that the separater isn't in the phenotype names
grep("_ALPHA_", parameters$phen, value = T)

# read in el data-o
alph_res <- lapply(parameters$file, function(file) {	
	# Read in results file
	dat <- read_ldak(paste0(path, file))
	return(dat)
})

# name the results
names(alph_res) <- unique_nam

save(alph_res, file = dat_name)

# load in the data
load(dat_name)

# extract the likelihood data
lik_res <- alph_res %>%
	map_df("likelihood_estimates", .id = "unique_id") %>%
	separate(col = unique_id, into = c("phen", "alpha"), sep = "_ALPHA_") %>%
	dplyr::filter(variable == "Alt_Likelihood")

# summarise likelihoods
lik_res %>%
	group_by(alpha) %>%
	summarise(med_lik = median(value))
anov <- aov(value ~ as.factor(alpha), data = lik_res)
summary(anov)
res.aov <- aov(weight ~ group, data = my_data)

cut_count <- lik_res %>% count(alpha)
phen_count <- lik_res %>% count(phen)

# Make some plots
p1_lik <- ggplot(lik_res, aes(x = alpha)) +
	geom_histogram(stat = "count") +
	labs(title = "Successful REML analyses at differing alpha values",
		x = "Alpha")

ggsave("results/sens/alpha_hist.pdf", plot = p1_lik)

# Extract the phenotypes where the analysis ran for all matrices
phen_out <- phen_count[phen_count$n == 9, "phen"][[1]]
res_out <- filter(lik_res, phen %in% phen_out)
dim(res_out)

summary(res_out$value)

# likelihood varying with alpha
p2_lik <- ggplot(lik_res, aes(x = alpha, y = value)) +
	geom_boxplot() +
	labs(x = "Alpha", 
		 y = "likelihood") +
	theme_bw()

ggsave("results/sens/alpha_likelihood_boxplot.pdf", plot = p2_lik)

########
#
# now with m2 estimates!
#
########

# extract m2 data
m2_res <- alph_res %>%
	map_df("her_estimates", .id = "unique_id") %>%
	separate(col = unique_id, into = c("phen", "alpha"), sep = "_ALPHA_") %>%
	dplyr::filter(Component == "Her_ALL") %>%
	rename(m2 = Heritability, se = Her_SD) %>%
	dplyr::select(phen, m2, se, alpha)


p1_m2 <- ggplot(m2_res, aes(x = as.factor(alpha), y = m2)) +
	geom_boxplot() +
	labs(title = "Association between alpha values and m2",
		x = "Alpha", 
		y = bquote(m^2)) +
	theme_bw()

ggsave("results/sens/alpha_m2_boxplot.pdf", plot = p1_m2)

temp_res <- m2_res %>%
	dplyr::filter(m2 < 0.3) %>%
	dplyr::filter(m2 > -0.3)

p1.5_m2 <- ggplot(temp_res, aes(x = as.factor(alpha), y = m2)) +
	geom_boxplot() +
	labs(title = "Association between alpha values and m2",
		x = "Alpha", 
		y = bquote(m^2))

ggsave("results/sens/alpha_m2_boxplot_zoomed.pdf", plot = p1.5_m2)

set.seed(2)
#phen_out <- phen_count[phen_count$nn == 8, "phen"]
res_out <- filter(m2_res, phen %in% sample(phen, 10))

p2_m2 <- ggplot(res_out, aes(x = as.factor(alpha), y = m2, group = phen, colour = phen)) +
	geom_line() +
	geom_point() +
	theme(legend.position = "none") +
	labs(title = "Association between alpha values and m2",
		x = "Alpha", 
		y = bquote(m^2))

ggsave("results/sens/alpha_m2_line.pdf", plot = p2_m2)

# Look at association of cutoff threshold and trend (assoc of m2 and EWAS hit number) 

# Read in some EWAS hits data at e-5
ewas_dat <- read.delim(paste0("../EWAS/results/", timepoints, "/derived/top_hits_e-5.txt"), stringsAsFactors = F)

# extract hits
phen_hits <- as.data.frame(table(ewas_dat$phen)) %>%
	rename(phen = Var1, hit_count = Freq) %>%
	dplyr::filter(phen %in% ind_vars)

# join hits to ldak res
hit_res <- left_join(phen_hits, m2_res) %>%
	mutate(p = pnorm(-abs(m2 / se))) %>%
	mutate(p = ifelse(m2 < 0, 1, p)) %>%
	mutate(var_exp = ifelse(p < 0.05, "Yes", "No")) %>%
	filter(!is.na(alpha)) 
head(hit_res)

p3_m2 <- ggplot(hit_res, aes(x = var_exp, y = log10(hit_count))) +
	geom_boxplot() +
	facet_grid(. ~ alpha) +
	labs(x = bquote(m^2 ~ "> 0"))

ggsave("results/sens/alpha_hit_count_box.pdf", plot = p3_m2)

head(arrange(hit_res, desc(hit_count)), n = 10)
head(arrange(hit_res[hit_res$var_exp == "Yes",], hit_count), n = 100)

# FIN