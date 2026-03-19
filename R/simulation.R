#' Simulate count data using fixed, parametric, or empirical gene effects
#'
#' @param n_genes number of genes to simulate
#' @param guides_per_gene number of guides per gene to simulate
#' @param n_control number of control samples
#' @param n_treatment number of treatment samples
#' @param p_effects proportion of genes with a non-zero effect
#' @param p_positive proportion of effects that are positive
#' @param fold_change absolute fold change for fixed mode; default mean for parametric mode
#' @param expression mean baseline expression level
#' @param guide_efficiency probability a guide has a non-zero effect
#' @param guide_sd standard deviation for guide-level effect noise
#' @param n_ntc number of non-targeting control guides
#' @param seed integer seed for reproducibility
#' @param gene_params named list for parametric mode (e.g. \code{list(mu = log(4), sd = 0.2)})
#' @param disp_params named list of dispersion trend coefficients (e.g. \code{list(a = 0.1, b = 4)})
#' @param effect_mode one of \code{"parametric"}, \code{"fixed"}, or \code{"empirical"}
#' @param counts optional matrix of raw counts; if provided, used to estimate baseline expression and dispersion
#' @param quantiles quantile bounds for selecting extreme effects in empirical mode (e.g. \code{c(0.05, 0.95)})
#' @return list with simulated counts, guide-to-gene mapping, sample design, and true parameter values
#' @export
make_playdough <- function(n_genes,
                           guides_per_gene,
                           n_control = 2,
                           n_treatment = 2,
                           p_effects = 0.1,
                           p_positive = 0.1,
                           fold_change = 4,
                           expression = 5,
                           guide_efficiency = 0.8,
                           guide_sd = 1,
                           n_ntc = 1000,
                           seed = 13,
                           gene_params = list(mu = log(fold_change), sd = 1),
                           disp_params = list(a = 0.1, b = 4),
                           effect_mode = c("parametric", "fixed", "empirical"),
                           counts = NULL,
                           quantiles = c(0.05, 0.95)) {
  set.seed(seed)
  effect_mode <- match.arg(effect_mode)

  n_samples <- n_control + n_treatment
  n_effects <- ceiling(p_effects * n_genes)
  n_guides <- n_genes * guides_per_gene + n_ntc

  # guide to gene mapping
  guide_ids <- seq_len(n_guides)
  target_gene_map <- rep(seq_len(n_genes), each = guides_per_gene)

  if (n_ntc > 0) {
    neg_gene_start <- n_genes + 1

    if (n_ntc >= guides_per_gene) {
      n_full_genes <- n_ntc %/% guides_per_gene
      n_leftover <- n_ntc %% guides_per_gene
      neg_gene_map <- rep(neg_gene_start:(neg_gene_start + n_full_genes - 1),
                          each = guides_per_gene)
    } else {
      n_full_genes <- 1
      n_leftover <- 0
      neg_gene_map <- rep(neg_gene_start, n_ntc)
    }

    if (n_leftover > 0) {
      neg_gene_map <- c(neg_gene_map, neg_gene_start:(neg_gene_start + n_leftover - 1))
    }

    gene_map <- c(target_gene_map, neg_gene_map)
  } else {
    gene_map <- target_gene_map
  }

  guide_to_gene <- data.frame(sgRNA = guide_ids, gene = gene_map)

  # baseline expression and dispersion
  if (!is.null(counts)) {
    if (effect_mode == "empirical") stop("Counts must be provided for empirical mode.")

    counts_mat <- as.matrix(counts)
    row_means <- pmax(rowMeans(counts_mat), 1e-6)
    row_vars <- matrixStats::rowVars(counts_mat)
    beta0_real <- log(row_means + 1)
    phi_real <- pmax((row_vars - row_means) / row_means^2, 1e-3)

    moment_indices <- sample(nrow(counts_mat), n_guides, replace = TRUE)
    beta0_g <- beta0_real[moment_indices]
    phi_g <- phi_real[moment_indices]
  } else {
    n_peak <- round(n_guides * 0.95)
    n_tail <- n_guides - n_peak
    combined <- sample(c(
      rnorm(n_peak, mean = expression, sd = 0.75),
      runif(n_tail, min = 0.1, max = expression)
    ))
    beta0_g <- pmin(abs(combined), 10)
    binding <- rbinom(n_guides, 1, p = guide_efficiency)
  }

  # gene-level effects
  gene_effect <- numeric(n_genes)
  n_neg <- round(n_effects * (1 - p_positive))
  signs <- sample(c(rep(-1, n_neg), rep(1, n_effects - n_neg)))

  if (effect_mode == "fixed") {
    gene_effect[seq_len(n_effects)] <- log(fold_change) * signs
    guide_sd <- rep(0, n_guides)

  } else if (effect_mode == "parametric") {
    gene_effect[seq_len(n_effects)] <- abs(rnorm(n_effects, gene_params$mu, gene_params$sd)) * signs
    guide_sd <- rlnorm(n_effects, log(guide_sd), 0.25)

  } else {
    lfc <- log(row_means + 1) - median(log(row_means) + 1)
    if (!is.null(quantiles)) {
      q <- quantile(lfc, probs = quantiles, na.rm = TRUE)
      lfc <- lfc[lfc <= q[1] | lfc >= q[2]]
    }
    gene_effect[seq_len(n_effects)] <- sample(lfc, n_effects, replace = TRUE)
  }

  # guide-level effects
  beta1_g <- numeric(n_guides)
  targeting_indices <- seq_len(n_genes * guides_per_gene)

  for (g in seq_len(n_effects)) {
    guides_idx <- which(gene_map == g)
    beta1_g[guides_idx] <- rnorm(length(guides_idx), mean = gene_effect[g], sd = guide_sd[g]) *
      binding[guides_idx]
  }

  # count matrix
  treatment_indices <- (n_control + 1):n_samples
  size_factors      <- runif(n_samples, 0.8, 1.2)

  mu_mat <- matrix(exp(beta0_g), nrow = n_guides, ncol = n_samples)
  mu_mat[targeting_indices, treatment_indices] <- exp(
    beta0_g[targeting_indices] + beta1_g[targeting_indices]
  )
  mu_mat <- mu_mat * matrix(size_factors, nrow = n_guides, ncol = n_samples, byrow = TRUE)

  # dispersion
  if (effect_mode != "empirical") {
    phi_g <- abs(disp_params$a + disp_params$b / rowMeans(mu_mat))
  }

  # simulate
  sim_counts <- matrix(
    rnbinom(n_guides * n_samples, mu = as.vector(mu_mat), size = as.vector(1 / phi_g)),
    nrow = n_guides,
    ncol = n_samples
  )

  rownames(sim_counts) <- guide_ids
  colnames(sim_counts) <- c(paste0("control", seq_len(n_control)),
                            paste0("treatment", seq_len(n_treatment)))

  sample_design <- data.frame(
    sample = colnames(sim_counts),
    design = as.factor(c(rep("control", n_control), rep("treatment", n_treatment)))
  )

  # take note of controls
  if (n_ntc > 0) {
    pseudo_genes <- (n_genes + 1):max(guide_to_gene$gene)
    controls <- data.frame(
      guides = guide_to_gene$sgRNA[guide_to_gene$gene %in% pseudo_genes],
      index = which(guide_to_gene$gene %in% pseudo_genes)
    )
  } else {
    controls <- NULL
  }

  list(
    data = list(
      counts   = sim_counts,
      row_data = guide_to_gene,
      col_data = sample_design,
      controls = controls
    ),
    truth = list(
      arguments = lapply(as.list(match.call())[-1], eval, envir = parent.frame()),
      significant_genes  = which(gene_effect != 0),
      significant_guides = which(gene_map %in% which(gene_effect != 0)),
      beta0 = beta0_g,
      beta1 = beta1_g,
      mu = gene_effect,
      phi = phi_g,
      tau = c(guide_sd, rep(0, n_genes - n_effects)),
      sf = size_factors,
      disp_params = disp_params
    )
  )
}
