# -------------------------------------------------------
# Making the covariate file
# -------------------------------------------------------

# Covars:
# - genetic PCs 1:10
# - cell counts
# - age in mothers/YoB

# Need to make separate files for discrete and quantitative covars 

# File should look like (no headers):
# FID IID covar1 covar2 covar3 ...
# --  --  ----    ----   ----
# --  --  ----    ----   ----
# --  --  ----    ----   ----
# --  --  ----    ----   ----

rm(list = ls())

# -------------------------------------------------------
# Setup
# -------------------------------------------------------
pkgs <- c("tidyverse", "stringr")
lapply(pkgs, require, character.only = TRUE)

devtools::load_all("~/repos/usefunc/")

#timepoints <- c("FOM", "antenatal")
#timepoints <- c("FOM")
#timepoints <- c("F7")
#timepoints <- c("cord")
#timepoints <- c("15up")
#timepoints <- "FOF"
args <- commandArgs(trailingOnly = TRUE)
samplesheet_file <- args[1]
cell_counts_file <- args[2]

timepoints <- as.character(args[1])
#phen <- as.character(args[2])

#phen <- "BMI"

# ARIES samplesheet
load(samplesheet_file)

# Cell counts
cell_counts <- read.table(cell_counts_file, header = T)

fam <- read.table(paste0("kinship/methyl_", paste(timepoints, collapse = "_"), ".fam"), header = F)

# -------------------------------------------------------
# Extract the covariates needed from the data
# -------------------------------------------------------

#Principal components
PC_dat <- read.table(paste0("pca/", timepoints, "_pcs.eigenvec"), sep = " ", header = F, stringsAsFactors = F) 
head(PC_dat)
colnames(PC_dat) <- c("FID", "IID", paste0(rep("PC", times = 20), 1:20))
PC_dat$ALN <- gsub("[A-Z]", "", PC_dat[["FID"]])
PC_dat <- dplyr::select(PC_dat, -IID, -FID)

colnames(samplesheet)[colnames(samplesheet) == "Sample_Name"] <- "IID"
colnames(fam) <- c("FID", "IID", "PID", "MID", "SEX", "PHENO")
fam$FID <- as.character(fam$FID)

PCs <- left_join(fam, PC_dat, by = c("FID" = "ALN"))
PCs <- dplyr::select(PCs, -PID, -MID, -SEX, -PHENO)
PCs <- PCs[, 1:12] # Keep PCs 1:10 
head(PCs)
nrow(PCs) - nrow(PCs[complete.cases(PCs), ]) # 108 people removed due to lack of genotype info

# Put the data together
dat <- Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "IID"), list(fam, samplesheet, cell_counts, PCs))
head(dat)
dim(dat)
dat <- dat[!is.na(dat$PC1), ]

# Sanity check
stopifnot(sum(dat$FID.x == dat$ALN) == nrow(dat))
stopifnot(sum(dat$FID.y == dat$ALN) == nrow(dat))

# Select the covariates for the analysis
q_covars <- c("Bcell", "CD4T", "CD8T", "Gran", "Mono", "NK", "age", paste0("PC", 1:10))

quant_covars <- dplyr::select(dat, ALN, IID, one_of(q_covars)) # Input quantitative covars
head(quant_covars)

write.table(quant_covars, paste0("phen/", timepoints, "/", paste(timepoints, collapse = "_"), ".qcovar"), quote = F, row.names = F, col.names = F)




