# ------------------------------------------------------
# Script to run EWAS using meffil --> Updated from Dr. Gemma Sharp's script
# ------------------------------------------------------


# pull in args
args	<- 	commandArgs(trailingOnly = TRUE)
samplesheet_file <- toString(args[1])
BorM    <-      toString(args[2]) #Betas or M-values (B or M)
TP <-  toString(args[3]) #Time point (cord or F7 or 15up or antenatal or FOM)
Covariates <-toString(args[4]) #list of covariates (eg: m_age,mum_uni,matsm,parity i.e. commas but no spaces or quotation marks)
WD<- toString(args[5]) #Working directory
split <- toString(args[6]) #Number to split by
# Testing script without args

Phenofile <- paste0("methyl_variance/phen/ALSPAC/data_", TP, ".txt")

if (Covariates == "none") {
	Covariates <- ""
}

print(Trait)
print(CellData)
print(CellAdj)
print(Phenofile)
print(BorM)
print(TP)
print(Covariates)
print(WD)
print(split)

#set working directory
setwd(WD)

pkgs <- c("foreign", "meffil", "tidyverse", "cowplot", "lme4", "GenABEL")
lapply(pkgs, require, character.only = T)

source("EWAS/R/new_ewas.r")

devtools::load_all("~/repos/usefunc")


#Load description of samples
load(samplesheet_file)
samplesheet <- subset(samplesheet, time_point == TP) %>%
	dplyr::filter(is.na(duplicate.rm))
if(TP != "antenatal" & TP != "FOM") {
	qletB <- samplesheet$ALN[which(samplesheet$QLET == "B")] #find alns for multiple pregnancies
	samplesheet <- samplesheet[-which(samplesheet$ALN %in% qletB), ] #remove multiple pregnancies
}
#load the methylation data
load(paste0("EWAS/", TP, "_filtered_450k_data_batch_removed.RData"))

#Load phenotype data (this should be stored in your working directory)
#Pheno<-read.dta(paste0(Phenofile,".dta"))
# phen_list <- read.delim(Trait, header = F, stringsAsFactors = F)
# phen_list <- as.character(phen_list[[1]])
# phen_list[grep("^\\d", phen_list)] <- paste0("X", phen_list[grep("^\\d", phen_list)])
# phen_list[grep("^_", phen_list)] <- paste0("X", phen_list[grep("^_", phen_list)])
# phen_list <- phen_list[!grepl("aln", phen_list)]
dat <- read.delim(Phenofile, stringsAsFactors = F)
phen_list <- colnames(dat)[!grepl("aln|qlet", colnames(dat))]
if (split == "NA") {
    split <- paste0("1,", length(phen_list))
}

split <- as.numeric(unlist(strsplit(split, ",")))
phen_list <- phen_list[split[1]:split[2]]

# Extract phenotypes from dat
Pheno <- dplyr::select(dat, aln, one_of(phen_list))

# Load covariates used in GCTA analysis
covars <- read.table(paste0("methyl_variance/phen/", TP, "/", TP, ".qcovar"))
colnames(covars) <- c("aln", "Sample_Name", "Bcell", "CD4T", "CD8T", "Gran", "Mono", "NK", "age", paste0("PC", 1:10))
head(covars)

#Add Sample_Name to Pheno (assuming Pheno contains aln)
Pheno <- merge(Pheno, samplesheet[, c("ALN", "Sample_Name")], by.x = "aln", by.y = "ALN") %>%
	left_join(covars) %>%
	dplyr::filter(!duplicated(Sample_Name))
# Convert to M-values if necessary
if(BorM == "M"){
	meth <- log2(meth / (1 - meth))
}

Covs <- colnames(covars)[-c(1,2)]
no_cov <- length(Covs)

# this loads many meffil components that are not present when just loading the package 
# with "library(meffil)". This is required to use the updated meffil.ewas function
# The updated function is to speed up the process
devtools::load_all("~/repos/meffil")

for (i in phen_list) {
	print(i)
    do_ewas(i, paste0("EWAS/results/", TP, "/"))
}
