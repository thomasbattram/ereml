# -------------------------------------------------------
# Using REML to assess the proportion of total variance across many traits explained by DNAm
# -------------------------------------------------------

# phen file (no headers!)
# FID IID pheno
# --  --  ---
# --  --  ---
# --  --  ---

# -------------------------------------------------------
# Setup
# -------------------------------------------------------

pkgs <- c("tidyverse")
lapply(pkgs, require, character.only = T)

devtools::load_all("~/repos/usefunc")

#timepoints <- c("FOM", "antenatal")
#timepoints <- c("FOM")
#timepoints <- c("F7")
#ldak_mod <- TRUE

args <- commandArgs(trailingOnly = TRUE)

timepoints <- as.character(args[1])
sens <- as.character(args[2])
ldak_mod <- as.logical(args[3])
print(ldak_mod)

setwd("methyl_variance")

kinship_files <- c(paste0("kinship/methyl_", paste(timepoints, collapse = "_")), paste0("kinship/PED_methyl_", paste(timepoints, collapse = "_")))
write.table(kinship_files, file = paste0("kinship/methyl_", paste(timepoints, collapse = "_"), "_kinships.txt"), quote = F, row.names = F, col.names = F, sep = "\n")

#sens <- 1
#sens <- 2
#sens <- FALSE 
if (sens == 1) {
	out <- seq(1.5, 3, 0.25)
	kinships <- paste0("sens/methyl_", timepoints, "_", out, "cutoff")
} else if (sens == 2) {
	out <- seq(-2, 0, 0.25)
	kinships <- paste0("sens/methyl_", timepoints, "_", out, "_alpha")
} else {
	kinships <- paste0("methyl_", paste(timepoints, collapse = "_"))
}

if (ldak_mod) kinships <- paste0("mod2/", kinships)

dat <- read_delim("phen/ALSPAC/data_FOM.txt", delim = "\t")

phens <- colnames(dat)
phens <- phens[!(grepl("aln", phens))]

var_var <- apply(dat, 2, var, na.rm = T)
no_var_rm <- names(var_var)[var_var == 0]
no_var_rm <- no_var_rm[!is.na(no_var_rm)] # still some NAs because there are binary variables here...
dat <- dat %>%
	dplyr::select(-one_of(no_var_rm))

phens <- phens[!phens %in% no_var_rm]

parameters <- expand.grid(phen = phens, kinship = kinships)
parameters$mgrm <- 1

# -------------------------------------------------------

dat <- dplyr::select(dat, aln, alnqlet, everything())

# Making parameters for LDAK
parameters$phenfile <- with(parameters, paste("phen/", timepoints, "/", phen, ".txt", sep=""))
parameters$kinshipfile <- with(parameters, paste("kinship/", kinship, sep=""))
parameters$outfile <- with(parameters, paste("estimates/", timepoints, "/", kinship, "_", phen, sep=""))
#parameters$covar <- with(parameters, paste("phen/", paste(timepoints, collapse = "_"), ".covar", sep = ""))
parameters$qcovar <- paste("phen/", timepoints, "/", paste(timepoints, collapse = "_", sep = ""), ".qcovar", sep = "") # May need to change...

head(parameters)

# Create phenofile function
createPhen <- function(parameters, row, dat, fam, pheno_col)
{
	p <- parameters[row, ]
	#if(file.exists(p$phenfile)) return(NULL)
	require(GenABEL)
	#d <- subset(alldat, phen==p$phen)
	d <- subset(dat, !duplicated(aln))
	value <- grep(paste0("^", p$phen, "$"), colnames(d), value = T)
	stopifnot(value == p$phen)
	f <- fam[, 1:2]
	d <- merge(d, f, by.x="aln", by.y="FID")
	d <- subset(d, select=c(aln, IID, get(value)))
	# Rank transform data
	if (!is.binary(d[[value]]) && !is.monomorphic(d[[value]])) {
		d[[value]] <- rntransform(d[[value]])
		print("Transformed")
	}
	d <- subset(d, !duplicated(aln))
	write.table(d, file=p$phenfile, row=F, col=F, qu=F)
	return(NULL)
}

# Read in the fam file 
fam <- read.table(paste0("kinship/methyl_", paste(timepoints, collapse = "_"), ".fam"), header = F)
colnames(fam) <- c("FID", "IID", "PID", "MID", "SEX", "PHENO")
head(fam)

# function to run REML using GCTA
runGcta <- function(parameters, row, mgrm=1)
{
	# Removed --covar because there are no discrete covariables now
	p <- parameters[row,]
	if(file.exists(paste(p$outfile, ".hsq", sep=""))) return(NULL)
	if (mgrm > 1) p$kinshipfile <- paste0(p$kinshipfile, "_kinships.txt")
	grm <- ifelse(mgrm>1, "mgrm", "grm")
	cmd <- paste("gcta --", grm, " ", p$kinshipfile, " --reml --reml-no-constrain --reml-no-lrt --pheno ", p$phenfile, " --qcovar ", p$qcovar, " --out ", p$outfile, " --thread-num 10 --reml-maxit 500", sep="")
	system(cmd)
}

# function to run REML using LDAK
runldak <- function(parameters, row, mgrm=1) {
	p <- parameters[row, ]
	p$outfile <- paste0("ldak_", p$outfile)
	# if(file.exists(paste(p$outfile, ".reml", sep=""))) return(NULL)
	if (mgrm > 1) p$kinshipfile <- paste0(p$kinshipfile, "_kinships.txt")
	grm <- ifelse(mgrm>1, "mgrm", "grm")
	cmd <- paste("kinship/ldak5.linux --grm ", p$kinshipfile, " --reml ", p$outfile, " --constrain NO --reml-iter 500 --kinship-details NO --pheno ", p$phenfile, " --covar ", p$qcovar)
	system(cmd)
}

parameters <- arrange(parameters, phen)
if (any(duplicated(parameters$outfile))) stop("DUPLICATED OUTFILES!!!!")

for (i in 1:nrow(parameters)) {
	print(i)
	print(parameters[i, "phen"])
	# if (file.exists(paste0("ldak_",parameters[i, "outfile"], ".reml"))) next
	createPhen(parameters, row = i, dat = dat, fam = fam, pheno_col = i + 2)
	# runGcta(parameters, row = i, mgrm = 1) # removed the GCTA section as no point in running both
	runldak(parameters, row = i, mgrm = 1)
	# stopifnot(file.exists(paste0("ldak_", parameters[i, "outfile"], ".reml")))
}

