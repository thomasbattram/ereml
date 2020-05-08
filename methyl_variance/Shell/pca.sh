# -------------------------------------------------------
# Generating PCs in ARIES
# -------------------------------------------------------

# Set directory
cd methyl_variance/pca

# Select timepoint
#timepoints="cord"
#timepoints="15up"
#timepoints="F7"
#timepoints="FOM"
timepoints=$1

# Extract all the ALNs from the "FOM" timepoint
grep ${timepoints} ../../ARIES_sample_ids.txt | awk ' { print $1$3" "$1$3 }' > ${timepoints}.txt

# fill these in!!
gendir=""

data=""

# Get snp list with no long range LD regions
awk -f ${gendir}/scripts/ld.awk ${data}.bim > nold.txt

# Get independent SNPs excluding any long range LD regions
plink --bfile $data --exclude nold.txt --indep 100 5 1.01 --out indep

# Calculate PCs
plink --bfile $data --keep ${timepoints}.txt --extract indep.prune.in --pca 20 --out ${timepoints}_pcs
