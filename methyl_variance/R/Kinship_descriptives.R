# -------------------------------------------------------
# Producing the PED-like GRM
# -------------------------------------------------------

pkgs <- c("tidyverse", "gridExtra")
lapply(pkgs, require, character.only = T)

devtools::load_all("~/repos/usefunc")

# source("R/GRM_functions.R")

args <- commandArgs(trailingOnly = TRUE)
samplesheet_file <- args[1] 
cell_counts_file <- args[2]
timepoints <- "FOM"

setwd("methyl_variance")

# Read in the methylation GRM
GRM_dat_b <- read_grm(rootname = paste0("kinship/methyl_", timepoints))
GRM_dat_g <- read_grm(rootname = paste0("kinship/mod2/methyl_", timepoints))
#GRM_dat <- readGRM(rootname = "/panfs/panasas01/sscm/tb13101/main_project/ALSPAC_EWAS/methyl_variance/kinship/methyl_FOM_antenatal")
#GRM_dat <- read_grm(rootname = "/panfs/panasas01/sscm/tb13101/main_project/ALSPAC_EWAS/methyl_variance/kinship/mod2/methyl_FOM")
grm_list <- list(blanket = GRM_dat_b, grouping = GRM_dat_g)
str(grm_list)

# -------------------------------------------------------
# Examine the data
# -------------------------------------------------------

extract_diag <- function(grm_dat) {
	stopifnot(is.list(grm_dat))
	diag <- filter(grm_dat$grm, id1 == id2) %>% 
		mutate(ALN = grm_dat$id$V1) %>%
		mutate(IID = grm_dat$id$V2)

	off_diag <- dplyr::filter(grm_dat$grm, id1 != id2)
	out_list <- list(diag = diag, off_diag = off_diag)
	return(out_list)
}

diag <- lapply(grm_list, extract_diag)

lapply(diag, function(x) {
	d <- var(x$diag$grm) 
	o <- var(x$off_diag$grm)
	return(list(diag=d, off_diag=o))
})

diag_grm_hist <- lapply(seq_along(diag), function(x) {
	diag_dat <- diag[[x]]$diag
	nam <- names(diag)[x]
	hist_p <- ggplot(diag_dat, aes(x = grm)) +
		geom_histogram(fill = "blue", colour = "black") +
		theme_bw() + 
		labs(title = nam)
	return(hist_p)
})

ggsave(paste0("kinship/descriptives/", timepoints, "_diag_hist.pdf"), 
	   plot = marrangeGrob(diag_grm_hist, nrow = 1, ncol = 1))

off_diag_grm_hist <- lapply(seq_along(diag), function(x) {
	diag_dat <- diag[[x]]$off_diag
	nam <- names(diag)[x]
	hist_p <- ggplot(diag_dat, aes(x = grm)) +
		geom_histogram(fill = "blue", colour = "black") +
		theme_bw() + 
		labs(title = nam)
	return(hist_p)
})

ggsave(paste0("kinship/descriptives/", timepoints, "_offdiag_hist.pdf"), 
	   plot = marrangeGrob(off_diag_grm_hist, nrow = 1, ncol = 1))


# looking at the max and minimum off diagonal values for blanket
# model
off_diag <- diag[["blanket"]]$off_diag

off_list <- list()
i=unique(off_diag$id1)[1]
for (i in unique(off_diag$id1)) {
	temp <- subset(off_diag, id1 == i | id2 == i)
	x <- extract_sum_stats(temp$grm)
	x$id <- i
	off_list[[i]] <- x
}


off_dat <- do.call(rbind, off_list)
head(off_dat)
g_off_dat <- off_dat %>%
	gather(key = "statistic", value = "value", -q1, -q3, -mean, -id)

p <- ggplot(g_off_dat, aes(x = value, fill = statistic)) +
	geom_histogram(position = "identity", alpha = 0.3)

ggsave("results/off_diagonal_max_and_min.pdf", plot = p)

# -------------------------------------------------------
# Examine outliers in the data
# -------------------------------------------------------

# Read in samplesheet
load(samplesheet_file)
samplesheet$aln <- samplesheet$ALN
samplesheet <- filter(samplesheet, time_point == "FOM")
# Cell counts
cell_counts <- read.table(cell_counts_file, header = T)

# Read in the phenotype data - smoking, BMI, alcohol consumption, 
smoke <- read.delim("../FOM_smoke_res.txt", stringsAsFactors = F)
BMI <- read.delim("phen/ALSPAC/BMI_Mum.txt", stringsAsFactors = F)
colnames(BMI)[3] <- "BMI"
alcohol <- read.delim("../FOM_alcohol_res.txt", stringsAsFactors = F)

diag$aln <- as.integer(diag$ALN)

diag_res <- Reduce(function(...) merge(..., all=TRUE, by="aln"), list(diag, smoke, BMI, alcohol, samplesheet))

diag_res <- left_join(diag_res, cell_counts, by = c("Sample_Name" = "IID"))
head(diag_res)
res <- dplyr::select(diag_res, aln, IID, grm, BMI, CIGS_smoked_per_day, Drinks_per_week, age, BCD_plate, Bcell, CD4T, CD8T, Gran, Mono, NK) %>%
	.[complete.cases(.), ]

# pdf(paste0("kinship/descriptives/", timepoints, "_diag_smoking_association.pdf"))
# plot(diag_res$grm, diag_res$CIGS_smoked_per_day)
# dev.off()

# What is the association of the diagonal values with the different terms
fit <- lm(grm ~ BMI + CIGS_smoked_per_day + Drinks_per_week + age + BCD_plate + Bcell + CD4T + CD8T + Gran + Mono + NK, data = res)
summary(fit)

res_hgrm <- filter(res, grm > 1.5)

fit2 <- lm(grm ~ BMI + CIGS_smoked_per_day + Drinks_per_week + age + BCD_plate + Bcell + CD4T + CD8T + Gran + Mono + NK, data = res_hgrm)
summary(fit2)

# Looking at each trait individually
summary(lm(grm ~ BMI, data = res_hgrm))
summary(lm(grm ~ CIGS_smoked_per_day, data = res_hgrm))
summary(lm(grm ~ Drinks_per_week, data = res_hgrm))
summary(lm(grm ~ age, data = res_hgrm))
summary(lm(grm ~ BCD_plate, data = res_hgrm))

# no large association there...

# -------------------------------------------------------
# Remove outliers in the data
# -------------------------------------------------------
# Use different cut-offs for outlier removal:
# - go from 1.5 - 3 in steps of 0.25 and add on with no outliers removed
# - write out tables that look like this:
# Timepoint  outlier_cutoff people_removed mean_diag median_diag var_diag
# ----        ---------      ---------      -------   ------      -------
# ----        ---------      ---------      -------   ------      -------
# ----        ---------      ---------      -------   ------      -------
# ----        ---------      ---------      -------   ------      -------
# ----        ---------      ---------      -------   ------      -------

outlier_vals <- seq(1.5, 3, by = 0.25)

diag_b <- diag$blanket

out_tab <- data.frame(
	timepoint = timepoints,
	n = nrow(GRM_dat_b$id),
	cutoff = NA,
	mean_diag = mean(diag_b$grm),
	median_diag = median(diag_b$grm),
	var_diag = var(diag_b$grm))

for (i in 2:(length(outlier_vals) + 1)) {
	out <- outlier_vals[i - 1]
	# Extract outliers
	outliers <- dplyr::filter(diag_b, grm >= out)

	new_GRM <- dplyr::filter(GRM_dat_b$grm, !(id1 %in% outliers$id1)) %>%
		dplyr::filter(!(id2 %in% outliers$id2))

	# Extract GRM diagonal
	diag_temp <- dplyr::filter(new_GRM, id1 == id2)

	# Write out the table
	out_tab[i, "timepoint"] <- timepoints
	out_tab[i, "n"] <- nrow(diag_temp)
	out_tab[i, "cutoff"] <- out
	out_tab[i, "mean_diag"] <- mean(diag_temp$grm)
	out_tab[i, "median_diag"] <- median(diag_temp$grm)
	out_tab[i, "var_diag"] <- var(diag_temp$grm)

	# Write out the new GRMs
	new_GRM_IDs <- dplyr::filter(GRM_dat_b$id, !(V2 %in% outliers$IID))

	GRM_out <- list(grm = new_GRM, id = new_GRM_IDs)

	write_grm(GRM_out, rootname = paste0("kinship/sens/methyl_", timepoints, "_", out, "cutoff"))

}

out_tab <- arrange(out_tab, desc(n))

grm_sum_stats_blanket <- extract_sum_stats(diag$blanket$diag$grm)
grm_sum_stats_grouping <- extract_sum_stats(diag$grouping$diag$grm)

out_list <- list(out_tab, grm_sum_stats_blanket, grm_sum_stats_grouping)
names(out_list) <- c("outlier_statistics_tab", "grm_diagonal_stats_blanket", "grm_sum_stats_grouping")

save(out_list, file = "results/grm_outlier_sum_stats.RData")

# write.table(out_tab, "results/grm_outlier_sum_stats.txt", quote = F, row.names = F, col.names = T, sep = "\t")

# Now take 10 random phenotypes (metabolites?) and test how these different values of cutoff influence the Vm value
# Can make a forest with Vm 95% CIs for each phenotype and each cutoff

# Removing diagonal elements above 2
non_outliers <- dplyr::filter(diag_b, grm < 2)
nrow(non_outliers) # 923 
outliers <- dplyr::filter(diag_b, grm >= 2)
nrow(outliers) # 17

var(non_outliers$grm) # 0.098


new_GRM <- dplyr::filter(GRM_dat_b$grm, !(id1 %in% outliers$id1)) %>%
	dplyr::filter(!(id2 %in% outliers$id2))
nrow(GRM_dat_b$grm) - nrow(new_GRM) 

new_GRM_IDs <- dplyr::filter(GRM_dat_b$id, !(V2 %in% outliers$IID))
nrow(new_GRM_IDs)

GRM_out <- list(grm = new_GRM, id = new_GRM_IDs)
str(GRM_out)


# Worth running REML with and without the outliers
write_grm(GRM_out, rootname = paste0("kinship/methyl_", timepoints, "_outliers_removed"))

