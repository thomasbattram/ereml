# -------------------------------------------------------
# Extracting ALSPAC data
# -------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
path <- args[1]
setwd(path)
# library(devtools)
# document()
# build()
# install()
#install.packages("devtools")
#library(devtools)
#install_github("thomasbattram/alspac")

pkgs <- c("alspac", "tidyverse", "haven", "readxl", "varhandle")
lapply(pkgs, require, character.only = T)
setDataDir("/Volumes/ALSPAC-Data")

devtools::load_all("~/Desktop/projects/Main_project/repos/usefunc")

data(current)
data(useful)
timepoints <- c("FOM")

# -------------------------------------------------------
# RUN THIS TO UPDATE THE DICTIONARIES
# -------------------------------------------------------
# current <- createDictionary("Current", name="current")
# useful <- createDictionary("Useful_data", name="useful")

# -------------------------------------------------------
# Filter out data not present enough in ARIES
# -------------------------------------------------------

# Read in the ARIES IDs and extract ones from timepoint of interest
IDs <- read.delim("ALSPAC_data/ARIES_sample_ids.txt", header = T, stringsAsFactors = F)
IDs <- dplyr::filter(IDs, time_point == timepoints)
str(IDs)

# # -------------------------------------------------------
# # For kinship descriptives
# # -------------------------------------------------------

# # Smoking variables for kinship descriptives stuff
# smok_res <- filt_res[, colnames(filt_res) %in% c("aln", "CIGS_smoked_per_day")]
# write.table(smok_res, file = "ALSPAC_EWAS/ALSPAC_data/FOM_smoke_res.txt", quote = F, row.names = F, col.names = T, sep = "\t")

# # smoking, alcohol use, BMI, cell count and batch variables 

# # Alcohol variables for kinship descriptives stuff
# alc_vars <- findVars("alcohol")

# filt_vars <- rbind(filter(alc_vars, cat3 %in% c("Adult", "Mother")), filter(alc_vars, cat4 %in% c("Adult", "Mother")))
# filt_vars <- filter(filt_vars, type %in% c("integer", "numeric"))
# tot_alc <- filter(filt_vars, lab == "Total Alcoholic Drinks Per Week")
# res <- extractVars(tot_alc)
# colnames(res)[2] <- "Drinks_per_week"
# filt_res <- filter(res, aln %in% IDs$ALN)
# dim(filt_res) # 896 people left

# write.table(filt_res, file = "ALSPAC_EWAS/ALSPAC_data/FOM_alcohol_res.txt", quote = F, row.names = F, col.names = T, sep = "\t")

# -------------------------------------------------------
# Extract all the Mum's data in ALSPAC
# -------------------------------------------------------

if (timepoints %in% c("FOM", "antenatel")) {
	tim <- c("Adult", "Mother")
} else if (timepoints %in% c("FOF")) {
	tim <- c("Adult", "Father")
} else {
	tim <- c("Child")
}


# Get adult data
unique(current$path)

# paths of interest 
PoI <- c(
	"Current/Other/Obstetric",
	"Current/Other/Samples/Mother",
	"Current/Other/Social_Class",
	"Current/Clinic/Adult",
	"Current/Quest/Mother", 
	"Current/Other/Social_Class"
)

new_current <- current %>%
	dplyr::filter(path %in% PoI)
sort(unique(new_current$obj))

# ------------------------------------------------------------------------------------
# Extract data 
# ------------------------------------------------------------------------------------

# extract variables!
result <- extractVars(new_current)

res <- result
# remove multi variables
mult_cols <- grep("mult", colnames(result), value = T)
grep("aln|qlet", colnames(res), value = TRUE)
# delete extra columns minus the aln, qlet and alnqlet columns
col_rm <- colnames(res)[!colnames(res) %in% new_current$name]
col_rm <- col_rm[!col_rm %in% c("qlet", "alnqlet")]
res <- res[,!colnames(res) %in% col_rm]
qlet_cols <- grep("qlet", colnames(res), value = T)
# remove the extra aln_qlet column and just keep FOM aries people
res <- res %>%
	dplyr::select(-aln_qlet) %>%
	dplyr::filter(aln %in% IDs$ALN)

# function to get all the descriptive names of the alspac variables
get_full_names <- function(input) {
	out <- map_chr(seq_along(input), function(x) {
		if (is.null(attributes(input[[x]]))) return(colnames(input[x]))
		return(attr(input[[x]], "label"))
	})
	return(out)
}

nams <- get_full_names(res)

all_labels <- map_df(seq_along(res), function(x) {
	if (is.null(attributes(res[[x]]))) return(NULL)
	labels <- attr(res[[x]], "labels")
	out <- data.frame(lab = names(labels), value = labels)
	return(out)
})

dim(res)

# ------------------------------------------------------------------------------------
# Start cleaning data
# ------------------------------------------------------------------------------------

grep("consent", all_labels$lab, value = TRUE, ignore.case = T)

fact_res <- as_factor(res)

cat_vars <- map_chr(seq_along(fact_res), function(x) {
	phen_levels <- levels(fact_res[[x]])
	if (is.null(phen_levels)) return("NULL")
	# remove missing and consent withdrawn
	phen_levels <- phen_levels[-grep("missing|consent", phen_levels, ignore.case = TRUE)]
	if (length(phen_levels) < 20 & length(phen_levels) > 4) {
		return(colnames(fact_res[x]))
	} else {
		return("NULL")
	}
})
cat_vars <- cat_vars[cat_vars != "NULL"]
cat_var_full_names <- get_full_names(res[, cat_vars])

# remove cat vars!!! 
fact_res2 <- fact_res %>%
	dplyr::select(-one_of(cat_vars))

# convert variables to missing
na_vars <- unique(c(grep("missing|consent", all_labels$lab, ignore.case = TRUE, value = TRUE),
				  "Not YE v5 / YHL v2", 
				  "YE short", 
				  "YHL", 
				  "HaB short", 
				  "FTG", 
				  "YP short", 
				  "NS/NK", 
				  "Not stated", 
				  "Not completed", 
				  "DNA", 
				  "NA", 
				  "NS/NA", 
				  "DK"))

fact_res2[] <- lapply(seq_along(fact_res2), function(x) {
	var <- fact_res2[[x]] 
	out <- mapvalues(var, from=na_vars, to=rep(NA, length(na_vars)))
	return(out)
})

# extract levels! 
few_levels <- map_df(seq_along(fact_res2), function(x) {
	if (is.null(attributes(fact_res2[[x]]))) return(NULL)
	labels <- attr(fact_res2[[x]], "levels")
	if (length(labels) > 10) return(NULL)
	out <- expand.grid(colnames(fact_res2[x]), labels)
	return(out)
})
dim(few_levels)

uniq_values <- few_levels %>%
	group_by(Var1) %>%
	summarise(n = n())

table(uniq_values$n)
uv4 <- uniq_values %>%
	dplyr::filter(n == 4)

few_levels %>%
	dplyr::filter(Var1 %in% uv4$Var1)

to_check <- c("a051", "j131c", "j132c", "j133c", "j134c", "j135c", "j136c", "j140a")

str(fact_res[, to_check])
levels(fact_res[, to_check[8]])
attributes(res[[to_check[2]]])
res[[to_check[2]]]
# a051 is definitely a categorical --> answers are yes, no, incoming calls
# others are uncertain. I.e. there are values in the variable that aren't labelled but the question
# should be a yes or no
# Removing all of them

to_rm <- uv4$Var1

uv3 <- uniq_values %>%
	dplyr::filter(n == 3)

few_levels %>%
	dplyr::filter(Var1 %in% uv3$Var1) %>%
	dplyr::filter(Var1 != "V4860A") # removing from view as labels are too long

to_rm <- c(to_rm, uv3$Var1, dplyr::filter(uniq_values, n==1)$Var1)

fin_fact_res <- fact_res2 %>%
	dplyr::select(-one_of(to_rm))

dim(fin_fact_res)

# --------------------------------------------------------------
# remove phenotypes with too much missing data
# --------------------------------------------------------------
missing_dat <- map_df(seq_along(fin_fact_res), function(x) {
	out <- data.frame(phen = colnames(fin_fact_res[x]), na_count = sum(is.na(fin_fact_res[[x]])))
	return(out)
})
sum(missing_dat$na_count > nrow(res)/2) # 698 phenotypes have over 50% missing data
to_rm <- missing_dat[missing_dat$na_count > nrow(res)/2, "phen"]
fin_fact_res <- fin_fact_res %>%
	dplyr::select(-one_of(to_rm))

# select only variables left
res2 <- fin_fact_res
res2[,-c(1:3)] <- lapply(seq_along(res2)[-c(1:3)], function(x) {
	print(x)
	var <- res2[[x]]
	label <- attributes(var)$label
	out <- unfactor(var)
	attributes(out)$label <- label
	return(out)
})
dim(res2)
# ------------------------------------------------------------------------------------
# save duplicated column names for later! 
# ------------------------------------------------------------------------------------
new_nams <- get_full_names(res2)
dup_labs <- new_nams[duplicated(new_nams)] # need to come back to these if we want repeated measures!
x <- which(new_nams %in% dup_labs)

res2 <- res2[, -x]
dim(res2)
# ------------------------------------------------------------------------------------
# Sorting binary vals
# ------------------------------------------------------------------------------------
res_bin <- res2[, colnames(res2) %in% uniq_values$Var1]
res2 <- res2[, !(colnames(res2) %in% colnames(res_bin))]

dim(res2)
res2[] <- lapply(seq_along(res2), function(x) {
	print(x)
	if (colnames(res[x]) %in% qlet_cols) return(res2[[x]])
	label <- attributes(res2[[x]])$label
	out <- as.numeric(res2[[x]])
	attributes(out)$label <- label
	return(out)
})

# Removal of categories where there are <100 values
missing <- sapply(res2, function(x) {sum(is.na(x))})
names(missing)
vars_rm <- missing[missing > (nrow(res2) / 2)]
length(vars_rm) # 212 variables removed due to lack of people
res3 <- dplyr::select(res2, -one_of(names(vars_rm)))
dim(res3)

uniq_vals <- sapply(res_bin, function(x) length(unique(x[!is.na(x)])))
# remove any without 2 unique values
to_rm <- names(uniq_vals)[uniq_vals != 2]
res_bin <- res_bin[, !(colnames(res_bin) %in% to_rm)]

# remove binary variables with too few cases
# -- removing if less than 10% cases or controls
few_cases <- apply(res_bin, 2, function(x) {sum(unique(x)[1] == x, na.rm = T) < (nrow(res_bin) / 10)})
few_controls <- apply(res_bin, 2, function(x) {sum(unique(x)[2] == x, na.rm = T) < (nrow(res_bin) / 10)})

cc_var_rm <- unique(c(names(which(few_controls)), names(which(few_cases))))
length(cc_var_rm) # 727
res_bin <- dplyr::select(res_bin, -one_of(cc_var_rm))

# remove binary variables with too much missing
missing <- sapply(res_bin, function(x) {sum(is.na(x))})
names(missing)
vars_rm <- missing[missing > (nrow(res_bin)/2)]
length(vars_rm) # 0 variables removed due to lack of people

res3 <- cbind(res3, res_bin)
dim(res3)
# ------------------------------------------------------------------------------------
# Finishing tidying data + saving it all
# ------------------------------------------------------------------------------------

# Swap the labels and the alspac names
res3_nam <- get_full_names(res3)
res4 <- map_dfc(seq_along(res3), function(x) {
	var <- res3[[x]]
	attributes(var)$alspac_name <- colnames(res3[x])
	return(var)
})
colnames(res4) <- res3_nam

# Rename the headings to remove all the unusable characters for GCTA
colnames(res4) <- gsub("\\%", "percent", colnames(res4))
colnames(res4) <- gsub("[[:punct:]]", "_", colnames(res4))
colnames(res4) <- gsub(" ", "_", trimws(colnames(res4)))

# CHECK FOR PHENOTYPES THAT ARE WHACK! 
colnames(res4)

# Total --> Unclear what this means, removing
# .*modes.*missing --> Unclear what this means, removing
# Hours_per_week --> Unclear what this means, removing
# .*modes_for_missing_values --> Unclear what this means, removing
# DV__Length_of_time_between_last_meal_and_blood_sample_taken__hours___FOM1 --> combined with one below! 
# DV__Length_of_time_between_last_meal_and_blood_sample_taken__minutes___FOM1
# .*fieldworker.* --> Removed, completely irrelevant except for with FOM1 blood-draw...
# DV__Length_of_time_between_last_meal_and_blood_sample_taken__hours___FOM2 --> Removed, completely irrelevant except for with FOM1 blood-draw...
# DV__Length_of_time_between_last_meal_and_blood_sample_taken__minutes___FOM2 --> Removed, completely irrelevant except for with FOM1 blood-draw...
# DV__Length_of_time_between_last_meal_and_blood_sample_taken__hours___FOM3 --> Removed, completely irrelevant except for with FOM1 blood-draw...
# DV__Length_of_time_between_last_meal_and_blood_sample_taken__minutes___FOM3 --> Removed, completely irrelevant except for with FOM1 blood-draw...
# Time_last_eaten__hour___FOM4 --> Removed, completely irrelevant except for with FOM1 blood-draw...
# Time_last_eaten__minutes___FOM4 --> Removed, completely irrelevant except for with FOM1 blood-draw...
# DV__Length_of_time_between_last_meal_and_blood_sample_taken__hours___FOM4 --> Removed, completely irrelevant except for with FOM1 blood-draw...
# DV__Length_of_time_between_last_meal_and_blood_sample_taken__minutes___FOM4 --> Removed, completely irrelevant except for with FOM1 blood-draw...
# DV__Maternal_bonding_score__modes_missing_values --> Unclear what this means, removing
# 1_gestage_days_ --> Unclear what this means, removing
# 1_gestage_weeks_ --> Unclear what this means, removing
# 2_```` --> Unclear what this means, removing
# 3_ --> Unclear what this means, removing
# . --> Unclear what this means, removing
# 14 --> Unclear what this means, removing
# Standard_error_of.* --> removing
# LCA_of_all_women__probability_of_membership_of_Class_1 --> Unclear what this means, removing
# DV__Gestation_in_days_based_on_LMP --> Unclear what this means, removing
# DV__Gestation_in_weeks_based_on_LMP --> Unclear what this means, removing
# Comments_made_on_B303 --> Unclear on who made the comments, removing
# Comments_made_on_B310 --> Unclear on who made the comments, removing
# Any_comments --> Unclear on who made the comments, removing
# Comments_on_birth --> Unclear on who made the comments, removing
# Form_number --> Unclear what this means, removing
# DV__Frequency_Mum_has_taken_vitamins_since_study_child_was_18_months_old__Y_N --> Keeping 
# Data_available__as_of_30_04_07_ --> Unclear what this means, removing
# Source__Single_babies_clean_files --> Unclear what this means, removing
# Data_available__as_of_31_12_07_ --> Unclear what this means, removing
# Data_available__as_of_18_05_12_ --> Unclear what this means, removing
# Questionnaire_completed --> Removing
# Questionnaire_data_available__as_of_30_04_05_ --> Unclear what this means, removing

# write out all the meta data!!! 
phen_list <- map_df(seq_along(res4), function(x) {
	if (colnames(res4[x]) %in% c(qlet_cols, "aln")) return(NULL)
	out <- data.frame(
		phen = colnames(res4[x]),
		binary = is.binary(res4[[x]]),
		n = sum(!is.na(res4[[x]])),
		alspac_name = attributes(res4[[x]])$alspac_name,
		unedited_label = attributes(res4[[x]])$label
		) %>%
		mutate(obj = new_current[new_current$name == alspac_name, "obj"])
	return(out)
})

# now check traits --> MANUALLY
res4$Total
phen_list[phen_list$phen == "Total",]
new_current[new_current$obj == "b_4f.dta","path"]
total_dat <- read_dta("/Volumes/ALSPAC-Data/Current/Quest/Mother/b_4f.dta")
total_dat$b915
###
grep(".*modes.*missing", colnames(res4), value = T)
res4$Affection_score___modes_missing_data
###
res4$Hours_per_week
phen_list[phen_list$phen == "Hours_per_week",]
new_current[new_current$obj == "e_4f.dta","path"]
hpw_dat <- read_dta("/Volumes/ALSPAC-Data/Current/Quest/Mother/e_4f.dta")
hpw_dat$e524b
colnames(hpw_dat)
hpw_dat$e524a ### this is hours per week the partner looks
###
grep("gestage_days_", colnames(res4), value = T)
res4[["1_gestage_days_"]]
phen_list[phen_list$phen == "1_gestage_days_",]
new_current[new_current$obj == "OA_r1b.dta","path"]
gest_dat <- read_dta("/Volumes/ALSPAC-Data/Current/Other/Obstetric/OA_r1b.dta")
###
res4$LCA_of_all_women__probability_of_membership_of_Class_1
###
res4$Comments_made_on_B303
###
res4$DV__Frequency_Mum_has_taken_vitamins_since_study_child_was_18_months_old__Y_N

#### remove variables that are getting removed
to_rm <- c("Total", grep(".*modes.*missing", colnames(res4), value = T), "Hours_per_week",
		   grep(".*fieldworker.*", colnames(res4), value = T), 
		   "DV__Length_of_time_between_last_meal_and_blood_sample_taken__hours___FOM2", 
		   "DV__Length_of_time_between_last_meal_and_blood_sample_taken__minutes___FOM2", 
		   "DV__Length_of_time_between_last_meal_and_blood_sample_taken__hours___FOM3", 
		   "DV__Length_of_time_between_last_meal_and_blood_sample_taken__minutes___FOM3", 
		   "Time_last_eaten__hour___FOM4", 
		   "Time_last_eaten__minutes___FOM4", 
		   "DV__Length_of_time_between_last_meal_and_blood_sample_taken__hours___FOM4", 
		   "DV__Length_of_time_between_last_meal_and_blood_sample_taken__minutes___FOM4", 
		   "DV__Maternal_bonding_score__modes_missing_values", 
		   grep("gestage_days_|gestage_weeks_", colnames(res4), value = T), 
		   grep("Standard_error_of.*", colnames(res4), value = T), 
		   "Comments_made_on_B303", 
		   "Comments_made_on_B310", 
		   "Any_comments", 
		   "Comments_on_birth", 
		   "Form_number", 
		   "Data_available__as_of_30_04_07_", 
		   "Source__Single_babies_clean_files", 
		   "Data_available__as_of_31_12_07_", 
		   "Data_available__as_of_18_05_12_", 
		   "Questionnaire_completed", 
		   "Questionnaire_data_available__as_of_30_04_05_")

res4 <- res4 %>%
	dplyr::select(-one_of(to_rm))

phen_list <- phen_list %>%
	dplyr::filter(!(phen %in% to_rm))

#### merge the blood sample data together:
# DV__Length_of_time_between_last_meal_and_blood_sample_taken__minutes___FOM1 
to_merge <- c("DV__Length_of_time_between_last_meal_and_blood_sample_taken__minutes___FOM1", 
			  "DV__Length_of_time_between_last_meal_and_blood_sample_taken__hours___FOM1")
res4 <- res4 %>%
	mutate(DV__Length_of_time_between_last_meal_and_blood_sample_taken_minutes_FOM1 = 
		DV__Length_of_time_between_last_meal_and_blood_sample_taken__minutes___FOM1 + 
		DV__Length_of_time_between_last_meal_and_blood_sample_taken__hours___FOM1 * 60) %>%
	dplyr::select(-one_of(to_merge))

df_temp <- data.frame(phen = "DV__Length_of_time_between_last_meal_and_blood_sample_taken_minutes_FOM1", 
					  binary = FALSE, 
					  n = phen_list[phen_list$phen == to_merge[1], "n"], 
					  alspac_name = NA, 
					  unedited_label = "DV__Length_of_time_between_last_meal_and_blood_sample_taken_minutes_FOM1", 
					  obj = NA)

phen_list <- phen_list %>%
	dplyr::filter(!phen %in% to_merge) %>%
	bind_rows(df_temp)

file_nam <- paste0("ALSPAC_data/phenotype_metadata_", timepoints, ".txt")
write.table(phen_list, file = file_nam, quote = F, col.names = T, row.names = F, sep = "\t")

file_nam <- paste0("ALSPAC_data/data_", timepoints, ".txt")
write.table(res4, file = file_nam, quote = F, col.names = T, row.names = F, sep = "\t")

# Set new password each time
PASSWORD <- ""

zip(gsub(".txt", ".zip", file_nam), 
    files = file_nam, 
    flags = paste("--password", PASSWORD))

system(paste0("rm ", file_nam))

# FIN! 
