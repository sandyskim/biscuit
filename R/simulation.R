#' simulate count data using either fixed, parametric, or empirical gene effects
#'
#' @param n_genes number of genes to simulate
#' @param guides_per_gene number of guides per gene to simulate
#' @param n_control number of control samples
#' @param n_treatment number of treatment samples
#' @param p_effects proportion of genes with a non-zero effect
#' @param p_positive proportion of effects that are positive
#' @param fold_change absolute fold change of non-zero effects for fixed-effect mode, default mean for parametric mode
#' @param gene_params optional list for parametric mode (e.g. list(mu = log(4), sd = 0.2))
#' @param guide_sd optional standard deviation for guide effect simulation
#' @param n_ntc number of non-targeting control guides
#' @param seed integer seed for reproducibility
#' @param effect_mode "fixed", "parametric", or "empirical": modes that determine how gene effects are simulated
#' @param counts optional matrix of raw experimental counts; if provided, will use to calculate baseline expression and guide-wise dispersion
#' @param quantiles optional quantiles (e.g. c(0.05, 0.95)) to draw significant effects from in empirical mode
#' @return list with simulated counts, guide to gene mapping, sample design and true parameter values
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
  n_samples  <- n_control + n_treatment
  n_effects  <- ceiling(p_effects * n_genes)
  n_guides <- (n_genes * guides_per_gene) + n_ntc

  # create guide to gene mapping
  guide_ids <- 1:n_guides
  target_gene_map <- rep(1:n_genes, each = guides_per_gene)

  # full genes
  n_full_genes <- n_ntc %/% guides_per_gene
  # leftover guides
  n_leftover <- n_ntc %% guides_per_gene

  neg_gene_start <- n_genes + 1

  # full pseudo-genes
  neg_gene_map <- rep(neg_gene_start:(neg_gene_start + n_full_genes - 1),
                      each = guides_per_gene)

  # distribute leftover guides across the first few pseudo-genes
  if (n_leftover > 0) {
    leftover_genes <- neg_gene_start:(neg_gene_start + n_leftover - 1)
    neg_gene_map <- c(neg_gene_map, leftover_genes)
  }

  # full gene map
  gene_map <- c(target_gene_map, neg_gene_map)

  # guide to gene mapping
  guide_to_gene <- data.frame(sgRNA = guide_ids, gene = gene_map)

  # simulate guide-wise dispersion and baseline expression
  if (!is.null(counts)) {
    if (effect_mode == 'empirical') {
      stop('no counts provided for empirical estimation')
    }

    # if provided with counts, calculate empirical estimates
    counts_mat <- as.matrix(counts)

    # method of moments
    row_means <- rowMeans(counts_mat)
    row_means <- pmax(row_means, 1e-6)
    row_vars  <- matrixStats::rowVars(counts_mat)

    beta0_real <- log(row_means + 1)
    phi_real   <- pmax((row_vars - row_means) / row_means^2, 1e-3)

    # randomly sample pairs of moments (beta0, phi)
    moment_indices <- sample(nrow(counts_mat), n_guides, replace = TRUE)
    beta0_g <- beta0_real[moment_indices]
    phi_g <- phi_real[moment_indices]
  }
  else {
    # beta0
    n_peak <- round(n_guides * 0.95)
    peak <- rnorm(n_peak, mean = expression, sd = 0.75)
    n_tail <- n_guides - n_peak
    tail <- runif(n_tail, min = 0.1, max = expression)
    combined <- sample(c(peak, tail))
    beta0_g <- abs(combined)
    binding <- rbinom(n_guides, 1, p = guide_efficiency)
  }

  # simulate gene-level effects
  gene_effect <- numeric(n_genes)

  if (effect_mode == "fixed") {
    # fixed gene effect (one numeric value)
    signs <- sample(c(rep(-1, round(n_effects * (1 - p_positive))),
                      rep(1, n_effects - round(n_effects * (1 - p_positive)))))

    gene_effect[1:n_effects] <- log(fold_change) * signs
    guide_sd <- rep(0, n_guides)

  } else if (effect_mode == "parametric") {
    # parametric gene effect (distribution with given parameters)
    gene_mu <- gene_params$mu
    gene_sd <- gene_params$sd

    signs <- sample(c(rep(-1, round(
      n_effects * (1 - p_positive))),
      rep(1, n_effects - round(
        n_effects * (1 - p_positive)))))

    gene_effect[1:n_effects] <- abs(rnorm(n_effects, mean = gene_mu, sd = gene_sd)) * signs
    guide_sd <- rlnorm(n_effects, log(guide_sd), 0.25)

  } else if (effect_mode == "empirical") {
    # empirical gene effect (distribution estimated from empirical data)
    lfc <- log(row_means + 1) - median(log(row_means) + 1) # empirical log-fold change, blind to sample conditions
    if (!is.null(quantiles)) {
      # take extreme effects based on defined quantiles
      q <- quantile(lfc, probs = quantiles, na.rm = TRUE)
      lfc <- lfc[(lfc <= q[1] | lfc >= q[2])]
    }
    gene_effect[1:n_effects] <- sample(lfc, n_effects, replace = TRUE)
    guide_sd <- guide_sd
  }

  # simulate guide-level effects
  beta1_g <- rep(0, times = n_guides) # non-targeting guides don't have an actual effect, only subject to technical noise
  for (g in 1:n_effects) {
    guides_idx <- which(gene_map == g)
    beta1_g[guides_idx] <- rnorm(
      length(guides_idx),
      mean = gene_effect[g],
      sd = guide_sd[g]
    ) * binding[guides_idx]
  }

  # calculate baseline mean counts
  mu_mat <- matrix(exp(beta0_g), nrow = n_guides, ncol = n_samples)

  # simulate size factors
  size_factors <- runif(n_samples, 0.8, 1.2)

  # apply treatment effect
  treatment_indices <- (n_control + 1):n_samples
  targeting_indices <- 1:(n_genes * guides_per_gene)

  mu_mat[targeting_indices, treatment_indices] <- exp(beta0_g[targeting_indices] + beta1_g[targeting_indices])

  mu_mat <- mu_mat * matrix(size_factors,
                            nrow = n_guides,
                            ncol = n_samples,
                            byrow = TRUE)

  if(effect_mode != 'empirical') {
    a <- disp_params$a
    b <- disp_params$b
    phi_g = abs(a + b/(rowMeans(mu_mat)))
  }

  #  calculate dispersion
  size_mat <- 1 / (phi_g)

  # negative binomial count simulation
  sim_counts <- matrix(
    rnbinom(
      n_guides * n_samples,
      mu   = as.vector(mu_mat),
      size = as.vector(size_mat)
    ),
    nrow = n_guides,
    ncol = n_samples
  )

  # prepare metadata
  rownames(sim_counts) <- guide_ids
  colnames(sim_counts) <- c(paste0("control", 1:n_control),
                            paste0("treatment", 1:n_treatment))

  sample_design <- data.frame(sample = colnames(sim_counts),
                              design = as.factor(c(rep("control", n_control), rep("treatment", n_treatment))))

  pseudo_genes <- (n_genes + 1):max(guide_to_gene$gene)
  controls <- data.frame(
    guides = guide_to_gene$sgRNA[guide_to_gene$gene %in% pseudo_genes],
    index  = which(guide_to_gene$gene %in% pseudo_genes)
  )

  playdough <- list(
    data = list(
      counts   = sim_counts,
      row_data = guide_to_gene,
      col_data = sample_design,
      controls = controls
    ),
    truth = list(
      arguments = lapply(as.list(match.call())[-1], eval, envir = parent.frame()),
      significant_genes = which(gene_effect != 0),
      significant_guides = which(gene_map %in% which(gene_effect != 0)),
      beta0  = beta0_g,
      beta1  = beta1_g,
      mu     = gene_effect,
      phi    = phi_g,
      tau = c(guide_sd, rep(0, n_genes-n_effects)),
      sf = size_factors,
      disp_params = if (exists("disp_params")) disp_params else NULL
    )
  )

  return(playdough)
}
