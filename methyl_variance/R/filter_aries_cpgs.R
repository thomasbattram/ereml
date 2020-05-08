# -------------------------------------------------------
# Filter ARIES betas
# -------------------------------------------------------

# This script performs QC on ARIES DNA methylation data and 
# extracts the samples wanted for the analysis.
# ARIES data version = v4

# -------------------------------------------------------
# Setup
# -------------------------------------------------------

pkgs <- c("tidyverse", "meffil", "lme4", "GenABEL")
sapply(pkgs, require, character.only = T)

devtools::load_all("~/repos/usefunc/")

args <- commandArgs(trailingOnly = TRUE)
samplesheet_file <- args[1]
betas_file <- args[2]
fdata_file <- args[3]
detp_file <- args[4]
timepoints <- "FOM"

# ARIES samplesheet
load(samplesheet_file) 
samplesheet <- subset(samplesheet, time_point == timepoints) %>%
		dplyr::filter(is.na(duplicate.rm))
if(timepoints != "antenatal" & timepoints != "FOM") {
	qletB <- samplesheet$ALN[which(samplesheet$QLET == "B")] #find alns for multiple pregnancies
	samplesheet <- samplesheet[-which(samplesheet$ALN %in% qletB), ] #remove multiple pregnancies
}

# ARIES betas
load(betas_file)
meth <- beta[, samplesheet$Sample_Name] #keep the samples that correspond to the time point you're interested in
rm(beta)

# fdata
load(fdata_file)

# detection p values
load(detp_file)
pvals <- detection.p[,samplesheet$Sample_Name] #keep the samples that correspond to the time point you're interested in
rm(detection.p)

# list of probes that aren't problematic according to zhou et al. 2017
zhou_list <- read.delim("EWAS/retain_from_zhou.txt", header = F)
zhou_list <- as.character(zhou_list[[1]])

# -------------------------------------------------------
# 
# -------------------------------------------------------

# Regress out batch from the methylation betas
filt_sample <- filter(samplesheet, Sample_Name %in% colnames(meth))
index <- match(colnames(meth), filt_sample$Sample_Name)
filt_sample <- filt_sample[index, ]

# Sanity check
stopifnot(all(filt_sample$Sample_Name == colnames(meth)))

batch_rm <- TRUE
if (batch_rm) {
	meth_resid <- apply(meth, 1, function(x) {resid(lmer(as.numeric(x) ~ (1 | filt_sample$BCD_plate)))})
	# rownames(meth_resid) <- filt_sample$Sample_Name
	meth <- t(meth_resid)
	colnames(meth) <- filt_sample$Sample_Name
	rm(meth_resid)
}

stopifnot(all(colnames(meth) == samplesheet[["Sample_Name"]]))

#load annotation data
annotation <- meffil.get.features("450k")

#Filter meth data (remove sex chromosomes and SNPs and probes with high detection P-values)
pvalue_over_0.05 <- pvals > 0.05
count_over_0.05 <- rowSums(sign(pvalue_over_0.05))
Probes_to_exclude_Pvalue <- rownames(pvals)[which(count_over_0.05 > ncol(pvals) * 0.05)]
XY <- as.character(annotation$name[which(annotation$chromosome %in% c("chrX", "chrY"))])
SNPs.and.controls <- as.character(annotation$name[-grep("cg|ch", annotation$name)])
annotation<- annotation[-which(annotation$name %in% c(XY, SNPs.and.controls, Probes_to_exclude_Pvalue)), ]
meth <- subset(meth, row.names(meth) %in% annotation$name)
paste("There are now ", nrow(meth), " probes")
paste(length(XY), "were removed because they were XY")
paste(length(SNPs.and.controls), "were removed because they were SNPs/controls")
paste(length(Probes_to_exclude_Pvalue), "were removed because they had a high detection P-value")
rm(XY, SNPs.and.controls, pvals, count_over_0.05, pvalue_over_0.05, Probes_to_exclude_Pvalue)

filtered_vars <- c("detection_p_values", "on_XY", "SNPs/controls")

if (batch_rm) filtered_vars <- c(filtered_vars, "batch_removed")

# Retain from zhou
meth <- meth[rownames(meth) %in% zhou_list, ]
dim(meth)

q <- rowQuantiles(meth, probs = c(0.25, 0.75), na.rm = T)
iqr <- q[, 2] - q[, 1]
too.hi <- which(meth > q[,2] + 3 * iqr, arr.ind=T)
too.lo <- which(meth < q[,1] - 3 * iqr, arr.ind=T)
if (nrow(too.hi) > 0) meth[too.hi] <- NA
if (nrow(too.lo) > 0) meth[too.lo] <- NA

samp <- colnames(meth)

# rank transform the methylation data...
meth <- apply(meth, 1, rntransform)
meth <- t(meth)
colnames(meth) <- samp

# double check samples being used... 
dim(meth)
num_na <- apply(meth, 2, function(x){sum(is.na(x))})
rem_samp <- which(num_na > (0.05 * nrow(meth)))
meth <- meth[, -rem_samp]
dim(meth)

print(paste0("Number of samples removed = ", length(rem_samp)))

if (batch_rm) {
	nam <- paste0("EWAS/", timepoints, "_filtered_450k_data_batch_removed.RData")
} else {
	nam <- paste0("EWAS/", timepoints, "_filtered_450k_data.RData")
}

save(meth, file = nam)

print("Saved filtered methylation data")

# load("~/main_project/ALSPAC_EWAS/EWAS/FOM_filtered_450k_data.RData")
# dim(meth)


