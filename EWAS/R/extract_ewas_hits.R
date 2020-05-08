# -------------------------------------------------------
# Extract top hits from EWAS
# -------------------------------------------------------

pkgs <- c("tidyverse", "GenABEL", "haven")
lapply(pkgs, require, character.only = T)

devtools::load_all("~/repos/usefunc")

timepoints <- "FOM"

setwd("EWAS")

zhou_list <- read.delim("retain_from_zhou.txt", header = F)
zhou_list <- as.character(zhou_list[[1]])

# path <- "EWAS/"
# phens <- read.delim(paste0(path, "EWAS_phens.txt"), stringsAsFactors = F, header = F)
# phens <- as.character(phens[[1]])
# phens[grep("^\\d", phens)] <- paste0("X", phens[grep("^\\d", phens)])
# phens[grep("^_", phens)] <- paste0("X", phens[grep("^_", phens)])

# phens
phen_dat <- read_tsv("../methyl_variance/phen/ALSPAC/data_FOM.txt")
phens <- colnames(phen_dat)[-grep("aln|qlet", colnames(phen_dat))]
non_norm <- read.delim("../methyl_variance/phen/ALSPAC/non_normal_traits.txt", header = F, stringsAsFactors=F)
non_norm <- non_norm[[1]]
phens <- phens[!phens %in% non_norm]
# -------------------------------------------------------
# Extract all hits at varying p value thresholds
# -------------------------------------------------------
extract_hits <- function(path, datafile, threshold) {
	ifelse (grepl(".RData", datafile), load(paste0(path, datafile)), return(NULL))
	name <- gsub(".RData", "", datafile)
	dat_list <- vector(mode = "list", length = length(threshold))
	for (i in 1:length(threshold)) {
		t <- threshold[i]
		dat_list[[i]] <- dplyr::filter(res, p.all < t) %>%
			dplyr::filter(probeID %in% zhou_list) %>%
			arrange(p.all) %>%
			mutate(phen = name) %>%
			mutate(p_threshold = t) %>%
			mutate(hit_count = nrow(.))
		print(paste0("Extracted ", nrow(dat_list[[i]]), " CpGs associated with ", name, " at p < ", t))
		if (nrow(dat_list[[i]]) == 0) {
			dat_list[[i]] <- data.frame(probeID = NA, coef.all = NA, se.all = NA, p.all = NA, phen = name, p_threshold = t, hit_count = 0)
		}	
	}
	dat <- do.call(rbind, dat_list)
	print(summary(dat))
	return(dat)
}

path <- paste0("results/", timepoints, "/")
files <- list.files(path = path)
files <- grep("RData", files, value = T)

new_files <- paste0(phens, ".RData")
sum(new_files %in% files)
files <- new_files
###################################### 
thresholds <- 10^seq(-7, -3.0, by = 0.5)
names(thresholds) <- paste0("e", log10(thresholds))

top_hit_list <- vector(mode = "list", length = length(files))
names(top_hit_list) <- files

hits_file <- paste0(path, "derived/phen_hits.txt")
ewas_res_file <- paste0(path, "derived/all_ewas_res.txt")
if (file.exists(hits_file)) {
	phen_hits <- read_delim(hits_file, delim = "\t") %>%
		dplyr::filter(phen %in% phens)
	ewas_res <- read_delim(ewas_res_file, delim = "\t") %>%
		dplyr::filter(phen %in% phens)
}

stopifnot(all(unique(ewas_res$phen) %in% unique(phen_hits$phen)))

top_hit_list <- lapply(files, function(x) {
	phen <- gsub(".RData", "", x)
	if (!phen %in% phen_hits$phen) {
		tryCatch(extract_hits(path, x, thresholds),
			error = function(e) {err_msg(e)})
	}
})

top_hits <- do.call(rbind, top_hit_list)
rownames(top_hits) <- NULL

# Coumt hits for each phenotype

x <- dplyr::select(top_hits, phen, p_threshold, hit_count)
hits <- x[!duplicated(x), ]

print("Hit count of 0 for this many variables: ")
print(sum(hits$hit_count == 0, na.rm=T))

if (exists("phen_hits")) {
	hits <- rbind(phen_hits, hits) %>%
		dplyr::filter(!is.na(phen))
}

# remove the non-normal shite! 
hits <- hits %>%
	dplyr::filter(!phen %in% non_norm)

write.table(hits, file = paste0(path, "derived/phen_hits.txt"), quote = F, col.names = T, row.names = F, sep = "\t")

# write out the actual data also! 
y <- dplyr::select(top_hits, -p_threshold, -hit_count)

if (exists("ewas_res")) {
	y <- rbind(ewas_res, y) %>%
		dplyr::filter(!is.na(phen))
}

# remove the non-normal shite! 
y <- y %>%
	dplyr::filter(!phen %in% non_norm)

write.table(y, file = paste0(path, "derived/all_ewas_res.txt"), quote = F, col.names = T, row.names = F, sep = "\t")



