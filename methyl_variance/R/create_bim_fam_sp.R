# -------------------------------------------------------
# Creating BIM, FAM and SP files for REML analysis using LDAK
# -------------------------------------------------------

# Files needed
# ARIES betas
# ARIES samplesheet
# fdata

args <- commandArgs(trailingOnly = TRUE)
samplesheet_file <- args[1]
fdata_file <- args[2]

#timepoints <- c("FOM", "antenatal")
#timepoints <- c("FOM")
#timepoints <- c("F7")
#timepoints <- "FOF"
timepoints <- as.character(args[1])

pkgs <- c("tidyverse", "lme4")
lapply(pkgs, require, character.only = T)

# ARIES betas
load(paste0("EWAS/", timepoints, "_filtered_450k_data_batch_removed.RData"))
beta <- meth
rm(meth)
stopifnot(ncol(beta) < nrow(beta))
# ARIES samplesheet
load(samplesheet_file)

# fdata
load(fdata_file)

# Need to remove individuals who have multiple read outs at the same time point

dat <- filter(samplesheet, time_point %in% timepoints)
fam <- dat %>%
	mutate(PID = 0) %>%
	mutate(MID = 0) %>%
	mutate(PHENO = 0) %>%
	dplyr::select(ALN, Sample_Name, PID, MID, Sex, PHENO) %>%
	dplyr::filter(Sample_Name %in% colnames(beta))

remove <- as.numeric(names(which(table(fam$ALN) != length(timepoints))))
print(paste0("Number of people removed = ", length(remove)))
fam <- dplyr::filter(fam, !(ALN %in% remove))

# Sanity check 
stopifnot(as.numeric(names(which(table(fam$ALN) != length(timepoints)))) == 0)

unique(fam$Sex)
fam[["Sex"]] <- ifelse(fam$Sex == "F", 2, 1)
head(fam)
print("fam file ready")

# FAM format
# FID IID PID MID SEX PHENO
# --  --  --  --  --  ---
# --  --  --  --  --  ---
# --  --  --  --  --  ---
# --  --  --  --  --  ---

bim <- fdata.new %>%
	mutate(gd = 0) %>%
	mutate(a1 = "A") %>%
	mutate(a2 = "T") %>%
	dplyr::select(CHR, TargetID, gd, COORDINATE_37, a1, a2) %>%
	dplyr::filter(TargetID %in% rownames(beta))

unique(bim[["CHR"]])
bim[["CHR"]] <- as.character(bim[["CHR"]])

# Pointless bit because X and Y CHR CpGs have been removed....
bim[bim[["CHR"]] == "X", "CHR"] <- "23"
bim[bim[["CHR"]] == "Y", "CHR"] <- "24"

bim[["CHR"]] <- as.numeric(bim[["CHR"]])
bim[["COORDINATE_37"]] <- as.numeric(bim[["COORDINATE_37"]])

# BIM format
# CHR SNP/CPG GD POS A1 A2
# --	--    -- --  -- --
# --	--    -- --  -- --
# --	--    -- --  -- --
# --	--    -- --  -- --

print("bim file ready")

createSp <- function(bim, fam, dat)
{
	stopifnot(nrow(bim) == nrow(dat))
	stopifnot(nrow(fam) == ncol(dat))
	bim2 <- bim[order(bim$CHR, bim$COORDINATE_37), ]
	index <- match(bim2$TargetID, bim$TargetID)
	dat2 <- dat[index, ]
	return(list(sp=dat2, bim=bim2, fam=fam))
}

dat <- as.data.frame(beta) %>%
	dplyr::select(one_of(fam[["Sample_Name"]]))

stopifnot(all(colnames(dat) == fam[["Sample_Name"]]))
colnames(dat) <- fam[["ALN"]]
stopifnot(all(colnames(dat) == fam[["ALN"]]))

# SP format 
#        SAMPLE1 SAMPLE2 SAMPLE 3
# CPG1
# CPG2
# CPG3
# CPG4
print("betas ready")

sp <- createSp(bim, fam, dat)
print("all 3 files created")


write.table(sp$sp, file = paste0("methyl_variance/kinship/methyl_", paste(timepoints, collapse = "_"), ".sp"), quote = F, col.names = F, row.names = F)
write.table(sp$bim, file = paste0("methyl_variance/kinship/methyl_", paste(timepoints, collapse = "_"), ".bim"), quote = F, col.names = F, row.names = F)
write.table(sp$fam, file = paste0("methyl_variance/kinship/methyl_", paste(timepoints, collapse = "_"), ".fam"), quote = F, col.names = F, row.names = F)

