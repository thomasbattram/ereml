# -------------------------------------------------------
# Post EWAS analysis
# -------------------------------------------------------
# -------------------------------------------------------
# load data
# -------------------------------------------------------

timepoints="FOM"
pkgs <- c("tidyverse", "cowplot", "gridExtra", "GenABEL", "lmtest", "haven", "boot",
		  "caret", "plotROC", "psychometric", "pROC", "RColorBrewer", "ggExtra", "ggpubr", 
		  "pscl")
lapply(pkgs, require, character.only = T)

devtools::load_all("~/repos/usefunc")

setwd("EWAS")

# Read in the initial m2 data
load(paste0("../methyl_variance/results/", timepoints, "_inf_mod_ldak_output.RData"))
inf_mod <- res
load(paste0("../methyl_variance/results/", timepoints, "_mod2_ldak_output.RData"))
mod2 <- res
rm(res)

zhou_list <- read.delim("retain_from_zhou.txt", header = F)
zhou_list <- as.character(zhou_list[[1]])

# Read in top hits at all p vals
ewas_res <- read_delim(paste0("results/", timepoints, "/derived/all_ewas_res.txt"), delim = "\t")
top_hits_e5 <- ewas_res %>%
	dplyr::filter(p.all < 1e-5)
top_hits_e7 <- ewas_res %>%
	dplyr::filter(p.all < 1e-7)

# write.table(top_hits_e7, "results/FOM/derived/top_hits_e7.txt", quote = F, col.names = T, row.names = F, sep = "\t")

# Read in hit counts for phens
phen_hits <- read.delim(paste0("results/", timepoints, "/derived/phen_hits.txt"), stringsAsFactors = F)
phen_hits$threshold <- phen_hits$p_threshold

# read in original data
dat <- read_tsv(paste0("~/main_project/ALSPAC_EWAS/methyl_variance/phen/ALSPAC/data_", timepoints, ".txt"))

# read in independent traits
load("../methyl_variance/pruned_traits.RData")
ind_vars <- rownames(temp_dat)
ind_vars <- ind_vars[ind_vars %in% ewas_res$phen]

# read in the filtered likelihood results
lik_res <- read_tsv("../methyl_variance/results/FOM_filtered_m2_likelihood_estimates.txt") %>%
	dplyr::filter(phen %in% ind_vars)

# -------------------------------------------------------
# Summary of EWAS findings
# -------------------------------------------------------
### for paper: 
get_ewas_sum_stats <- function(dat, p) {
	ewas_sum <- dat %>%
		dplyr::filter(threshold == p)
	out <- list()
	
	out$n_cpg <- sum(ewas_sum$hit_count)
	out$n_trait <- sum(ewas_sum$hit_count > 0) 
	out$trait <- dplyr::filter(ewas_sum, hit_count == max(hit_count))

	return(out)

}

ind_phen_hits <- phen_hits %>%
	dplyr::filter(phen %in% ind_vars)

test <- get_ewas_sum_stats(phen_hits, 1e-7)
test2 <- get_ewas_sum_stats(ind_phen_hits, 1e-7) # for paper
test3 <- get_ewas_sum_stats(ind_phen_hits, 1e-7 / length(ind_vars)) # for paper
test4 <- get_ewas_sum_stats(ind_phen_hits, 1e-5) # for paper

# =======
sum_ewas_findings <- function(dat, threshold) {
	ewas_sum_res <- list()
	sum_assoc <- dat %>%
		group_by(phen) %>%
		summarise(assoc_num = length(probeID))

	ewas_sum_res$phen_num <- nrow(sum_assoc)
	ewas_sum_res$total_n <- sum(sum_assoc$assoc_num)
	
}

# for the paper:
ewas_sum_res <- list(p5 = list(), p7 = list())
# -- number of traits associated with >0 CpGs at p<10-7 and p<10-5
ewas_sum_res$p5$traits <- unique(top_hits_e5$phen)
ewas_sum_res$p7$traits <- unique(top_hits_e7$phen)
# -- what is the average number of associations
sum_assoc <- top_hits_e7 %>%
	group_by(phen) %>%
	summarise(assoc_num = length(probeID))

# Do traits that have a high Vm tend to have a higher number of hits associated with them?
lik_res_p <- lik_res %>%
	dplyr::filter(variable == "LRT_P") %>%
	dplyr::select(-variable) %>%
	dplyr::rename(p = value)

ldak_res <- list(blanket = inf_mod %>%
		map_df("her_estimates", .id = "phen"), 
	grouping = mod2 %>%
		map_df("her_estimates", .id = "phen")) %>%
	bind_rows(.id = "model") %>% 
	dplyr::filter(Component == "Her_ALL") %>%
	dplyr::select(model, phen, Heritability, Her_SD) %>%
	dplyr::rename(m2 = Heritability, se = Her_SD) %>%
	left_join(lik_res_p) %>%
	dplyr::filter(se > 0) %>%
	dplyr::filter(!is.na(p)) # removing NAs because that means alternative likelihood was lower than the null likelihood

# add fdr
ldak_res[ldak_res$model == "blanket", "fdr"] <- p.adjust(ldak_res[ldak_res$model == "blanket", "p", drop = T], method = "fdr")
ldak_res[ldak_res$model == "grouping", "fdr"] <- p.adjust(ldak_res[ldak_res$model == "grouping", "p", drop = T], method = "fdr")

ldak_res %>%
	group_by(model) %>%
	summarise(med_m2 = median(m2), mean_m2 = mean(m2))

ldak_res %>%
	dplyr::filter(fdr < 0.15)

ldak_res %>%
	dplyr::filter(fdr < 0.05)

# summary(ldak_res)

hit_res_all <- left_join(phen_hits, ldak_res) %>%
	dplyr::filter(!is.na(m2)) %>%
	dplyr::filter(phen %in% ind_vars)

##
# 
# for papier
#
##

hit_p5 <- hit_res_all %>%
	dplyr::filter(p_threshold == 1e-5) %>%
	mutate(low_m2 = ifelse(m2 < 0, "low", "normal"))

mycolours <- c("low" = "grey50", "normal" = "black")

mods <- unique(hit_p5$model)

ep1 <- lapply(mods, function(mod) {
	hit_dat <- hit_p5 %>%
		dplyr::filter(model == mod)
	p <- ggplot(hit_dat, aes(x = m2, y = hit_count)) +
  		geom_point(aes(colour = low_m2)) + 
  		scale_color_manual("Status", values = mycolours) +
  		labs(x = bquote(m^2), y = "Number of DMPs") +
  		geom_vline(xintercept = 0, linetype = "dashed", colour = "grey") +
  		theme(panel.background = element_rect(fill = "white"), 
    		  axis.line = element_line(size = 0.5, linetype = "solid", colour = "black"),
    		  legend.position = "none", axis.text = element_text(size = 10), 
    		  axis.title = element_text(size = 12))
  	pdf("temp.pdf")
  	p <- ggMarginal(p, type = "histogram", xparams = list(bins = 50), yparams = list(bins = 50))
  	dev.off()
  	return(p)
})
system("rm temp.pdf")

pdf("plots/m2_hit_count_scatter_p5.pdf")
print(marrangeGrob(ep1, ncol = 1, nrow = 2, top=NULL))
dev.off()

### testing!!! 

p <- ggplot(hit_p5, aes(x = m2, y = hit_count)) +
  		geom_point(aes(colour = low_m2)) + 
  		scale_color_manual("Status", values = mycolours) +
  		labs(x = bquote(m^2), y = "Number of DMPs") +
  		geom_vline(xintercept = 0, linetype = "dashed", colour = "grey") +
  		facet_grid(model~.) +
  		theme(panel.background = element_rect(fill = "white"), 
    		  axis.line = element_line(size = 0.5, linetype = "solid", colour = "black"),
    		  legend.position = "none", axis.text = element_text(size = 12), 
    		  axis.title = element_text(size = 14), strip.text = element_text(size = 14))
  		

ggsave("plots/m2_hit_count_scatter_p5_test.pdf", plot = p)

# -------------------------------------------------------
# Association between number of DMPs and m2 --> Searching for the best model! 
# -------------------------------------------------------
phen_dat <- dplyr::select(dat, one_of(hit_res_all$phen))

hit_res <- hit_res_all[hit_res_all$threshold == 1e-05, ] %>%
	dplyr::filter(model == "grouping") %>%
	mutate(var_exp = p < 0.05)

p <- ggplot(hit_res, aes(x = var_exp, y = log10(hit_count + 1), fill = model)) +
	geom_boxplot()

ggsave(file = paste0("plots/", timepoints, "hit_count_box.pdf"), plot = p)

head(arrange(hit_res, desc(hit_count)), n = 10)
head(arrange(hit_res[hit_res$var_exp == T,], desc(hit_count)), n = 100)

hit_res <- hit_res %>%
	dplyr::filter(Vm > -1)

t_hit_count <- log(hit_res$hit_count, base = 100)
hist(t_hit_count)

p2 <- ggplot(hit_res, aes(x = hit_count)) +
	geom_histogram(binwidth = 1) +
	theme_classic()
ggsave("plots/FOM_hit_count_distribution.pdf", plot = p2)

# hit count is too right skewed with too many 0 values to transform

# Linear regression assumptions
fit <- with(hit_res, lm(hit_count ~ m2 - 1))
summary(fit)
x <- resid(fit)
hist(x)
plot(x)
xy <- x[x < 10 & x > -10]
length(x) - length(xy) # 99
hist(xy)

# Poisson regression might be better! 

### testing the different Poisson regression methods:
# -- Ordinary Poisson
# -- Quasi-Poisson --> can estimate the dispersion (variability) in the data - fixed at one for ordinary Poisson
# -- Negative binomial --> more formal way to accommodate over-dispersion
# Two models that address additional zeros in the count data:
# -- Hurdle
# -- Zero-inflated Poisson

####### Ordinary Poisson
pois <- glm(hit_count ~ m2, family = poisson, data = hit_res)
summary(pois)

new_dat <- hit_res %>%
	dplyr::select(m2, var_exp)

# We test for goodness-of-fit of the model with a chi-square test based on the residual deviance
# and degrees of freedom.

1 - with(summary(pois), pchisq(deviance, df.residual))
1 - pchisq(summary(pois)$deviance, summary(pois)$df.residual)

# P = 0, therefore it is not a good fit!! 
dat <- cbind(new_dat,
	  Mean = predict(pois, newdata = new_dat, type = "response"), 
      SE = predict(pois, newdata = new_dat, type = "response", se.fit = T)$se.fit
      )

dat %>%
	group_by(var_exp) %>%
	summarise(pred_mean = median(Mean), pred_se = median(SE))

####### Quasi-Poisson
qpois <- glm(hit_count ~ m2, family = quasipoisson, data = hit_res)
summary(qpois)
# dispersion parameter taken to be 21.86!!! - very high!!
1 - pchisq(summary(qpois)$deviance, 
           summary(qpois)$df.residual
           )
# P = 0

####### Negative binomial
nbin <- glm.nb(hit_count ~ m2, data = hit_res)
summary(nbin)
1 - pchisq(summary(nbin)$deviance,
           summary(nbin)$df.residual
           )
# P = 0.121 --> much better

####### Hurdle
hurdlenb <- hurdle(hit_count ~ m2, data = hit_res, dist = "negbin")
summary(hurdlenb)
sum(predict(hurdlenb, type = "prob")[,1]) # by design it predicts number of 0s exactly!

hurdle_ori <- hurdle(hit_count ~ m2, data = hit_res)
summary(hurdle_ori)
AIC(hurdlenb)
AIC(hurdle_ori)

####### zero inflation negative binomial regression
zinb <- zeroinfl(hit_count ~ m2, data = hit_res, dist = "negbin")
summary(zinb)

####### zero inflation Poisson regression
zipr <- zeroinfl(hit_count ~ m2, data = hit_res)
summary(zipr)


fm <- list("ML-Pois" = pois, "Quasi-Pois" = qpois, "NB" = nbin, 
			"Hurdle-NB" = hurdlenb, "Hurdle" = hurdle_ori, "ZINB" = zinb, "ZIP" = zipr)
sapply(fm, function(x) coef(x))

model_res <- rbind(logLik = sapply(fm, function(x) round(logLik(x), digits = 0)),
	  Df = sapply(fm, function(x) attr(logLik(x), "df")))

model_res <- t(model_res) %>%
	as.data.frame() %>%
	rownames_to_column(var = "model")

write.table(model_res, "m2_hits_model_testing.txt", quote = F, row.names = F, col.names = T, sep = "\t")

## hurdle-nb seems to fit the data best! 
coef(hurdlenb)
confint(hurdlenb)
out <- cbind(coef(hurdlenb), confint(hurdlenb))
out <- as.data.frame(out) %>%
	rownames_to_column(var = "model") # for paper! 
colnames(out)[2] <- "estimate"

write.table(out, "dmp_m2_association_stats.txt", quote = F, row.names = F, col.names = T, sep = "\t")

# -------------------------------------------------------
# Ability of m2 to predict number of top hits
# -------------------------------------------------------

# Get number of CpGs tested
path <- "results/FOM/"
files <- list.files(path)
load(paste0(path, files[1]))
dplyr::filter(res, probeID %in% zhou_list) %>%
	nrow(.) -> ncpg
rm(res)

# define expected hit number
thresholds <- 10^(seq(-7, -3, by = 0.5))
exp_hits <- thresholds * ncpg
exp_h <- data.frame(threshold = thresholds, exp_hit_count = exp_hits)

for (i in thresholds) {
	hit_res_all$threshold <- ifelse(near(i, hit_res_all$threshold), i, hit_res_all$threshold)
}

roc_dat <- left_join(hit_res_all, exp_h)
head(roc_dat)
summary(hit_res_all$hit_count)
length(unique(hit_res_all$phen))

# make new variables saying whether the expected number of hits has been surpased or not
roc_dat$hits_above_exp <- with(roc_dat, hit_count > exp_hit_count)

##############


cols <- brewer.pal(n = length(unique(roc_dat$threshold)), name = "Set1")

test <- roc_dat %>%
	group_by(threshold) %>%
	summarise(n=n(), sig=sum(hits_above_exp))

test

pdf("brewer_test.pdf")
display.brewer.all()
dev.off()
# use code below 

roc_list_gcta <- list()
roc_list_ldak <- list()
auc_dat <- matrix(ncol = 5)
i=1
for (i in 1:length(unique(roc_dat$threshold))) {
	t <- unique(roc_dat$threshold)[i]
	logt <- as.character(-log10(t))
	temp_dat_gcta <- dplyr::filter(roc_dat, threshold == t) %>%
		dplyr::filter(model == "blanket")
	out <- with(temp_dat_gcta, roc(hits_above_exp ~ m2))
	roc_list_gcta[[logt]] <- out
	# auc <- gsub(".*\\:", "", out$auc)	
	auc <- cbind(t(as.data.frame(ci.auc(out))), logt)
	auc <- cbind(auc, "gcta")

	temp_dat_ldak <- dplyr::filter(roc_dat, threshold == t) %>%
		dplyr::filter(model == "grouping")
	out <- with(temp_dat_ldak, roc(hits_above_exp ~ m2))
	roc_list_ldak[[logt]] <- out
	# auc <- gsub(".*\\:", "", out$auc)
	# auc_dat[i, "ldak"] <- as.numeric(auc)
	auc2 <- cbind(t(as.data.frame(ci.auc(out))), logt)
	auc2 <- cbind(auc2, "ldak")
	auc <- rbind(auc, auc2)
	auc_dat <- rbind(auc_dat, auc)
}
rownames(auc_dat) <- NULL
auc_dat <- as.data.frame(auc_dat[-1,])
colnames(auc_dat) <- c("ci_low", "estimate", "ci_upper", "-log10(threshold)", "model")

auc_dat[, -c(4,5)] <- apply(auc_dat[, -c(4,5)], 2, function(x) {make_pretty(as.numeric(x), 2)})
auc_dat <- as.data.frame(auc_dat)

# auc plot
p_auc <- ggplot(auc_dat, aes(x = `-log10(threshold)`, y = estimate, fill = model)) +
	geom_bar(stat = "identity", position = position_dodge()) +
	geom_errorbar(aes(ymin = ci_low, ymax = ci_upper), position = position_dodge(0.9), width = 0.2)

ggplot2::ggsave("plots/auc_comparison.pdf", plot = p_auc)

auc_plot <- ggplot(auc_dat, aes(x = `-log10(threshold)`, y = estimate, colour = model)) +
  geom_point(position = position_dodge(width = 0.9)) +
  geom_linerange(aes(ymin = ci_low, ymax = ci_upper), position = position_dodge(width = 0.9)) +
  labs(x = bquote("P value threshold for EWAS (-log"[10] ~ "(P))"), 
       y = "Area under the curve")

auc_plot <- auc_plot +
    scale_colour_discrete(name="Model",
                         breaks=c("gcta", "ldak"),
                         labels=c("Blanket", "Grouping"))

ggplot2::ggsave("auc_plot_test3.pdf", plot = auc_plot)

# save data for paper
save(roc_list_gcta, roc_list_ldak, auc_dat, file = "model_comparison_roc_data.RData")

# make roc curve graph
### removing p = 1e-7 and p = 1e-6.5 as the number of hits above expected is much lower than the others
new_roc_list_gcta <- roc_list_gcta[-c(1,2)]
new_roc_list_ldak <- roc_list_ldak[-c(1,2)]
new_auc_dat <- auc_dat %>%
	dplyr::filter(!`-log10(threshold)` %in% c(6.5, 7))


p1 <- pROC::ggroc(new_roc_list_gcta) +
	geom_abline(intercept = 1, slope = 1, colour = "black", alpha = 0.6) +
	labs(title = "Blanket") +
	scale_colour_manual(values = cols) +
	theme_bw()

p2 <- pROC::ggroc(new_roc_list_ldak) +
	geom_abline(intercept = 1, slope = 1, colour = "black", alpha = 0.6) +
	labs(title = "Grouping") +
	scale_colour_manual(values = cols) + 
	theme_bw() + 
	theme(legend.position = "none", axis.text = element_text(size = 10))

# cols <- brewer.pal()

leg <- cowplot::get_legend(p1 + scale_colour_manual(values = cols, name = "-log10(P)\n[Blanket AUC, Grouping AUC]",
                         breaks = names(new_roc_list_gcta),
                         labels = paste0(names(new_roc_list_ldak), " [", dplyr::filter(new_auc_dat, model == "gcta")$estimate, ", ", dplyr::filter(new_auc_dat, model == "ldak")$estimate, "]")) +
						 		theme(legend.text = element_text(size = 8), legend.title = element_text(size = 8))
		)
p1 <- p1 + theme(legend.position = "none", axis.text = element_text(size = 10))


ptab <- ggtexttable(auc_dat[,-4], rows = NULL)
glist <- lapply(list(p1,p2), ggplotGrob)
class(ptab)
ggsave("plots/roc_plot.pdf", plot = marrangeGrob(glist, right = leg, nrow = 2, ncol = 1, top = NULL))
# for paper
print("ROC curves done")
