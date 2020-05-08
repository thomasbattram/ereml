# -------------------------------------------------------
# Extracting data from ldak
# -------------------------------------------------------

# -------------------------------------------------------
# Setup
# -------------------------------------------------------
setwd("methyl_variance")

pkgs <- c("tidyverse")
lapply(pkgs, require, character.only = T)

devtools::load_all("~/repos/usefunc")

# timepoints <- "FOM"
args <- commandArgs(trailingOnly = T)
timepoints <- args[1]

# Read in the phenotypes
path <- "phen/ALSPAC/"
#phens <- ifelse("FOM" %in% timepoints || "antenatal" %in% timepoints, read.delim(paste0(path, "phen_list", "_Mum.txt"), stringsAsFactors = F), read.delim(paste0(path, "phen_list", "_Child.txt"), stringsAsFactors = F))
dat <- read_tsv(paste0(path, "data_FOM.txt"))
phens <- colnames(dat)[-grep("aln|qlet", colnames(dat))]

# non normal phenotypes
non_norm <- read.delim("methyl_variance/phen/ALSPAC/non_normal_traits.txt", header = F, stringsAsFactors=F)
non_norm <- non_norm[[1]]
phens <- phens[!phens %in% non_norm]

# -------------------------------------------------------
# Extract data
# -------------------------------------------------------

mods <- c("inf_mod", "mod2")

mod = "inf_mod"
j = 1

# Independent vars only!!!
load("pruned_traits.RData")
ind_vars <- rownames(temp_dat)

# long loop to read in the data and condense it 
lapply(mods, function(mod) {
	# change path depending on model
	if (mod == "inf_mod") {
		path <- paste0("ldak_estimates/", timepoints, "/")	
	} else {
		path <- paste0("ldak_estimates/", timepoints, "/mod2/")
	}
	# get the files
	files <- list.files(path = path)
	files <- files[grep("methyl_", files, useBytes = TRUE)]
	dat_name <- paste0("results/", timepoints, "_", mod, "_ldak_output.RData")

	phen_files <- paste0("methyl_", timepoints, "_" , phens, ".reml")
	phen_files <- phen_files[phen_files %in% files]
	phens <- gsub("methyl_FOM_", "", phen_files)
	phens <- gsub(".reml", "", phens)
	# phens <- phens[phens %in% ind_vars]
	# loop over all the phenotypes --> read them in and keep as a long list
	res <- lapply(phens, function(phen) {
		file <- paste0("methyl_", timepoints, "_" , phen, ".reml")
		
		dat <- read_ldak(paste0(path, file))
		return(dat)
	})
	names(res) <- phens

	save(res, file = dat_name)
	return(NULL)
})

# -------------------------------------------------------
# filter out traits that have alternative likelihood lower than the null likelihood
# -------------------------------------------------------

load(paste0("results/FOM_inf_mod_ldak_output.RData"))
inf_mod <- res
load(paste0("results/FOM_mod2_ldak_output.RData"))
mod2 <- res
rm(res)

lik_res <- list(blanket = inf_mod %>%
		map_df("likelihood_estimates", .id = "phen"), 
	grouping = mod2 %>%
		map_df("likelihood_estimates", .id = "phen")) %>%
	bind_rows(.id = "model")

mod="blanket"
filtered_lik_res <- lapply(unique(lik_res$model), function(mod) {
	df <- lik_res %>%
		dplyr::filter(model == mod) %>%
		dplyr::filter(variable %in% c("Null_Likelihood", "Alt_Likelihood")) %>%
		spread(variable, value) %>%
		mutate(alt_minus_null = Alt_Likelihood - Null_Likelihood) %>%
		dplyr::filter(sign(alt_minus_null) != -1)

	out <- lik_res %>%
		dplyr::filter(model == mod) %>%
		dplyr::filter(phen %in% df$phen)
	return(out)
})

lik_res <- bind_rows(filtered_lik_res)

write.table(lik_res, file = "results/FOM_filtered_m2_likelihood_estimates.txt", 
			quote = F, row.names = F, col.names = T, sep = "\t")


