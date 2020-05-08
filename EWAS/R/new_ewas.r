# Adapted from meffil

meffil.ewas_new <- function (beta, variable, covariates = NULL, batch = NULL, weights = NULL, 
    cell.counts = NULL, isva = T, sva = T, n.sv = NULL, isva0 = F, 
    isva1 = F, winsorize.pct = 0.05, outlier.iqr.factor = NA, 
    most.variable = min(nrow(beta), 50000), featureset = NA, 
    random.seed = 20161123, verbose = F, crude = FALSE) {
    if (isva0 || isva1) 
        stop("isva0 and isva1 are deprecated and superceded by isva and sva")
    if (is.na(featureset)) 
        featureset <- guess.featureset(rownames(beta))
    features <- meffil.get.features(featureset)
    stopifnot(length(rownames(beta)) > 0 && all(rownames(beta) %in% 
        features$name))
    stopifnot(ncol(beta) == length(variable))
    stopifnot(is.null(covariates) || is.data.frame(covariates) && 
        nrow(covariates) == ncol(beta))
    stopifnot(is.null(batch) || length(batch) == ncol(beta))
    stopifnot(is.null(weights) || is.numeric(weights) && (is.matrix(weights) && 
        nrow(weights) == nrow(beta) && ncol(weights) == ncol(beta) || 
        is.vector(weights) && length(weights) == nrow(beta) || 
        is.vector(weights) && length(weights) == ncol(beta)))
    stopifnot(most.variable > 1 && most.variable <= nrow(beta))
    stopifnot(!is.numeric(winsorize.pct) || winsorize.pct > 0 && 
        winsorize.pct < 0.5)
    original.variable <- variable
    original.covariates <- covariates
    if (is.character(variable)) 
        variable <- as.factor(variable)
    stopifnot(!is.factor(variable) || is.ordered(variable) || 
        length(levels(variable)) == 2)
    msg("Simplifying any categorical variables.", verbose = verbose)
    variable <- simplify.variable(variable)
    if (!is.null(covariates)) 
        covariates <- do.call(cbind, lapply(covariates, simplify.variable))
    sample.idx <- which(!is.na(variable))
    if (!is.null(covariates)) 
        sample.idx <- intersect(sample.idx, which(apply(!is.na(covariates), 
            1, all)))
    msg("Removing", ncol(beta) - length(sample.idx), "missing case(s).", 
       verbose = verbose)
    if (is.matrix(weights)) 
        weights <- weights[, sample.idx]
    if (is.vector(weights) && length(weights) == ncol(beta)) 
        weights <- weights[sample.idx]
    beta <- beta[, sample.idx]
    variable <- variable[sample.idx]
    if (!is.null(covariates)) 
        covariates <- covariates[sample.idx, , drop = F]
    if (!is.null(batch)) 
        batch <- batch[sample.idx]
    if (!is.null(cell.counts)) 
        cell.counts <- cell.counts[sample.idx]
    if (!is.null(covariates)) {
        pos.var.idx <- which(apply(covariates, 2, var, na.rm = T) > 
            0)
        msg("Removing", ncol(covariates) - length(pos.var.idx), 
            "covariates with no variance.", verbose = verbose)
        covariates <- covariates[, pos.var.idx, drop = F]
    }
    if (crude == FALSE) {
    	covariate.sets <- list()
    } else {
    	covariate.sets <- list(none = NULL)
    }
    #covariate.sets <- list(none = NULL)
    if (!is.null(covariates)) 
        covariate.sets$all <- covariates
    if (is.numeric(winsorize.pct)) {
        msg(winsorize.pct, "- winsorizing the beta matrix.", 
            verbose = verbose)
        beta <- winsorize(beta, pct = winsorize.pct)
    }
    too.hi <- too.lo <- NULL
    if (is.numeric(outlier.iqr.factor)) {
        q <- rowQuantiles(beta, probs = c(0.25, 0.75), na.rm = T)
        iqr <- q[, 2] - q[, 1]
        too.hi <- which(beta > q[, 2] + outlier.iqr.factor * 
            iqr, arr.ind = T)
        too.lo <- which(beta < q[, 1] - outlier.iqr.factor * 
            iqr, arr.ind = T)
        if (nrow(too.hi) > 0) 
            beta[too.hi] <- NA
        if (nrow(too.lo) > 0) 
            beta[too.lo] <- NA
    }
    if (isva || sva) {
        beta.sva <- beta
        autosomal.sites <- meffil.get.autosomal.sites(featureset)
        autosomal.sites <- intersect(autosomal.sites, rownames(beta.sva))
        if (length(autosomal.sites) < most.variable) {
            warning("Probes from the sex chromosomes will be used to calculate surrogate variables.")
        }
        else {
            beta.sva <- beta.sva[autosomal.sites, ]
        }
        var.idx <- order(rowVars(beta.sva, na.rm = T), decreasing = T)[1:most.variable]
        beta.sva <- impute.matrix(beta.sva[var.idx, , drop = F])
        if (!is.null(covariates)) {
            cov.frame <- model.frame(~., data.frame(covariates, 
                stringsAsFactors = F), na.action = na.pass)
            mod0 <- model.matrix(~., cov.frame)
        }
        else mod0 <- matrix(1, ncol = 1, nrow = length(variable))
        mod <- cbind(mod0, variable)
        if (isva) {
            msg("ISVA.", verbose = verbose)
            set.seed(random.seed)
            isva.ret <- isva(beta.sva, mod, ncomp = n.sv, verbose = verbose)
            if (!is.null(covariates)) 
                covariate.sets$isva <- data.frame(covariates, 
                  isva.ret$isv, stringsAsFactors = F)
            else covariate.sets$isva <- as.data.frame(isva.ret$isv)
            cat("\n")
        }
        if (sva) {
            msg("SVA.", verbose = verbose)
            set.seed(random.seed)
            sva.ret <- sva(beta.sva, mod = mod, mod0 = mod0, 
                n.sv = n.sv)
            if (!is.null(covariates)) 
                covariate.sets$sva <- data.frame(covariates, 
                  sva.ret$sv, stringsAsFactors = F)
            else covariate.sets$sva <- as.data.frame(sva.ret$sv)
            cat("\n")
        }
    }
    analyses <- sapply(names(covariate.sets), function(name) {
        msg("EWAS for covariate set", name, verbose = verbose)
        covariates <- covariate.sets[[name]]
        ewas(variable, beta = beta, covariates = covariates, 
            batch = batch, weights = weights, cell.counts = cell.counts, 
            winsorize.pct = winsorize.pct)
    }, simplify = F)
    p.values <- sapply(analyses, function(analysis) analysis$table$p.value)
    coefficients <- sapply(analyses, function(analysis) analysis$table$coefficient)
    rownames(p.values) <- rownames(coefficients) <- rownames(analyses[[1]]$table)
    for (name in names(analyses)) {
        idx <- match(rownames(analyses[[name]]$table), features$name)
        analyses[[name]]$table$chromosome <- features$chromosome[idx]
        analyses[[name]]$table$position <- features$position[idx]
    }
    list(class = "ewas", version = packageVersion("meffil"), 
        samples = sample.idx, variable = original.variable[sample.idx], 
        covariates = original.covariates[sample.idx, , drop = F], 
        winsorize.pct = winsorize.pct, outlier.iqr.factor = outlier.iqr.factor, 
        most.variable = most.variable, p.value = p.values, coefficient = coefficients, 
        analyses = analyses, random.seed = random.seed, too.hi = too.hi, 
        too.lo = too.lo)
}

do_ewas <- function(phen, path) {
    if (file.exists(paste0(path, phen, ".RData"))) return(NULL)
    print(phen)
    
    #Prepare phenotype data
    temp_phen <- na.omit(Pheno[, c("Sample_Name", phen, Covs)])

    #Match meth to Pheno
    temp_meth <- meth[, na.omit(match(temp_phen$Sample_Name, colnames(meth)))]
    temp_phen <- temp_phen[match(colnames(temp_meth), temp_phen$Sample_Name), ]

    if (!(is.binary(temp_phen[, colnames(temp_phen) == phen]) || is.monomorphic(temp_phen[, colnames(temp_phen) == phen]))) {
            temp_phen[, colnames(temp_phen) == phen] <- rntransform(temp_phen[, colnames(temp_phen) == phen])
            print("Trait transformed")
    }

    ifelse(all(temp_phen$Sample_Name == colnames(temp_meth)), "meth and phenotype data successfully matched :) ", "Data not matched :(")
    temp_phen <- droplevels(temp_phen) # get rid of any empty factor levels 

    paste("There are ", nrow(temp_phen), " people in this analysis")
    "Here's a summary of the phenotype data:"
    summary(temp_phen)

    #Run EWAS using meffil
    tryCatch({
        obj <- meffil.ewas_new(temp_meth, variable = temp_phen[, 2], covariates = temp_phen[, -c(1,2)], winsorize.pct = NA, isva = F, sva = F, most.variable = min(nrow(temp_meth), 20000), outlier.iqr.factor = 3, verbose = TRUE, crude = FALSE)

        res <- data.frame(
        probeID = rownames(obj$analyses$all$table),
        coef.all = obj$analyses$all$table$coefficient,
        se.all = obj$analyses$all$table$coefficient.se,
        p.all = obj$analyses$all$table$p.value
        )

        save(res, file = paste0(path, phen, ".RData"))
        print(paste0("Results for ", phen, " saved."))
    }, error = function(e) {print(paste0("Error in EWAS of ", i, ". Variance of ", i, " = ", var(temp_phen[, 2])))})

}
