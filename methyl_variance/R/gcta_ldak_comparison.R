# --------------------------------------------------
# comparing blanket and grouping model
# --------------------------------------------------

setwd("methyl_variance")

pkgs <- c("tidyverse", "GenABEL", "RColorBrewer", "gridExtra", "ggExtra", "ggrepel")
lapply(pkgs, require, character.only = T)

devtools::load_all("~/repos/usefunc")

timepoints <- "FOM"

# load the data
load(paste0("results/", timepoints, "_inf_mod_ldak_output.RData"))
inf_mod <- res
load(paste0("results/", timepoints, "_mod2_ldak_output.RData"))
mod2 <- res
rm(res)

load("pruned_traits.RData")
ind_vars <- colnames(temp_dat)

# filtered likelihoods 
lik_res <- read_tsv("results/FOM_filtered_m2_likelihood_estimates.txt") %>%
	dplyr::filter(phen %in% ind_vars)

# --------------------------------------------------
# Checking the weights when using the grouping model
# --------------------------------------------------

weightings <- read_delim("kinship/mod2/methyl_FOM/weights.all", delim = " ")
sum(weightings$Weight)
sum(weightings$Weight != 0) # for paper
summary(weightings) # for paper
# --------------------------------------------------
# Testing differences between the likelihoods 
# --------------------------------------------------

lik_res2 <- lik_res %>%
	dplyr::filter(variable == "Alt_Likelihood") %>%
	spread(model, value) %>%
	mutate(blanket_div_grouping = blanket / grouping)

# check the differences between the likelihoods
wil_res <- wilcox.test(lik_res2$blanket, lik_res2$grouping, paired = TRUE)
wil_res # for paper. 
median(lik_res2$grouping) - median(lik_res2$blanket)

dif <- lik_res2$grouping - lik_res2$blanket
median(dif)

# shapiro.test(lik_res2$blanket)
# t_res <- t.test(lik_res2$blanket, lik_res2$grouping, paired = TRUE)
# t_res

table(sign(lik_res2$blanket - lik_res2$grouping))

write.table(lik_res2, file = "results/FOM_reml_likelihood_output.txt", quote = F, row.names = F, col.names = T, sep = "\t")

p <- ggplot(lik_res2, aes(x = blanket, y = grouping)) +
	geom_point() +
	geom_abline(colour = "red")

ggsave("results/plots/mod_comparison_likelihood_scatter.pdf", plot = p)

# check the differences between the LRT statistics
lik_res3 <- lik_res %>%
	dplyr::filter(variable == "LRT_Stat") %>%
	spread(model, value) %>%
	mutate(blanket_div_grouping = blanket / grouping)

wil_res <- wilcox.test(lik_res3$blanket, lik_res3$grouping, paired = TRUE)
wil_res # for paper. 
median(lik_res3$grouping) - median(lik_res3$blanket)

dif <- lik_res3$grouping - lik_res3$blanket
median(dif)


# check differences between null likelihoods -- should be none!
lik_res4 <- lik_res %>%
	dplyr::filter(variable == "Null_Likelihood") %>%
	spread(model, value) %>%
	mutate(blanket_div_grouping = blanket / grouping)

summary(lik_res4$blanket_div_grouping)

lik_res5 <- lik_res %>%
	dplyr::filter(variable == "LRT_P") %>%
	dplyr::select(-variable) %>%
	dplyr::rename(p = value)

# --------------------------------------------------
# Testing differences between the m2 values
# --------------------------------------------------

her_res <- list(blanket = inf_mod %>%
		map_df("her_estimates", .id = "phen"), 
	grouping = mod2 %>%
		map_df("her_estimates", .id = "phen")) %>%
	bind_rows(.id = "model") %>%
	dplyr::filter(Component == "Her_ALL") %>%
	mutate(m2 = Heritability) %>%
	mutate(se = Her_SD) %>%
	dplyr::select(model, phen, m2, se) %>%
	left_join(lik_res5) %>%
	dplyr::filter(se > 0) %>%
	dplyr::filter(!is.na(p)) # removing NAs because that means alternative likelihood was lower than the null likelihood

#########
p1 <- ggplot(her_res, aes(x = model, y = m2)) +
	geom_boxplot()

p2 <- ggplot(her_res, aes(x = model, y = -log10(se))) +
	geom_boxplot()

p3 <- ggplot(her_res, aes(x = model, y = -log10(p))) +
	geom_boxplot()

ggsave("results/plots/mod_comparison_box.pdf", plot = marrangeGrob(list(p1,p2,p3), nrow = 1, ncol = 3))

p4 <- ggplot(her_res, aes(x = m2, fill = model)) +
	geom_histogram(alpha = 0.2, position = "identity", bins = 100) +
	theme(legend.position = "none")

p5 <- ggplot(her_res, aes(x = -log10(se), fill = model)) +
	geom_histogram(alpha = 0.2, position = "identity", bins = 100) +
	theme(legend.position = "none")

p6 <- ggplot(her_res, aes(x = -log10(p), fill = model)) +
	geom_histogram(alpha = 0.2, position = "identity", bins = 100) +
	theme(legend.position = "none")

ggsave("results/plots/mod_comparison_hist.pdf", plot = marrangeGrob(list(p4,p5,p6), nrow = 3, ncol = 1))

#### for paper!!! --> find better way to do it!
blanket_fdr_res <- her_res %>%
	dplyr::filter(model == "blanket") %>%
	mutate(fdr = p.adjust(p, method = "BH")) %>%
	arrange(fdr)
dplyr::filter(blanket_fdr_res, fdr < 0.05)

grouping_fdr_res <- her_res %>%
	dplyr::filter(model == "grouping") %>%
	mutate(fdr = p.adjust(p, method = "BH")) %>%
	arrange(fdr)
dplyr::filter(grouping_fdr_res, fdr < 0.05)

#

her_res2 <- her_res %>%
	dplyr::select(model, phen, m2) %>%
	spread(model, m2) %>%
	dplyr::filter(!is.na(blanket)) %>%
	dplyr::filter(!is.na(grouping)) %>%
	mutate(highlight = ifelse(phen == "Smoked_CIGS_REG", "highlight", "normal"))

her_res2[her_res2$phen == "Smoked_CIGS_REG", "phen"] <- "Smoked_cigs_reg"

mycolours <- c("highlight" = "red", "normal" = "black")

p <- ggplot(her_res2, aes(x = blanket, y = grouping, colour = highlight)) +
	geom_point(aes(colour = highlight)) + 
	scale_color_manual("Status", values = mycolours) +
	geom_text_repel(data = her_res2[her_res2$phen == "Smoked_cigs_reg", ], aes(label = phen)) +
	labs(x = bquote(m^2 ~ " (blanket)"), y = bquote(m^2 ~ " (grouping)")) +
	geom_abline(intercept = 0, slope = 1, colour = "blue", linetype = "dashed") + 
	theme_bw()
p <- p + theme(legend.position = "none", axis.text = element_text(size = 12), 
			   axis.title = element_text(size = 14))

ggsave("results/plots/model_m2_comparison.pdf", plot = ggMarginal(p, type = "histogram", xparams = list(bins = 50), yparams = list(bins = 50)))


########## for paper (whole of bit below)
summary(her_res2)

diff_res <- her_res2 %>%
	mutate(dif = blanket - grouping) %>%
	arrange(desc(abs(dif)))

cor(diff_res$blanket, diff_res$grouping) # for paper

t.test(diff_res$blanket, diff_res$grouping, paired = T) # for paper
mean(diff_res$grouping) - mean(diff_res$blanket) # for paper
table(sign(diff_res$blanket - diff_res$grouping))

diff_res[2,"phen"]

her_res3 <- her_res %>%
	mutate(lower_ci = m2 - (1.96 * se)) %>%
	mutate(upper_ci = m2 + (1.96 * se)) 

her_res3 %>%
	dplyr::filter(phen == "Ratio_of_18_2_linoleic_acid_to_total_fatty_acids__percent___FOM2")

