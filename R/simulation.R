#' simulate count data using either fixed, parametric, or empirical gene effects
#'
#' @param n_genes number of genes to simulate
#' @param guides_per_gene number of guides per gene to simulate
#' @param n_control number of control samples
#' @param n_treatment number of treatment samples
#' @param p_effects proportion of genes with a non-zero effect
#' @param p_positive proportion of effects that are positive
#' @param fold_change absolute fold change of non-zero effects for fixed-effect mode, default mean for parametric mode
#' @param n_ntc number of non-targeting control guides
#' @param seed integer seed for reproducibility
#' @param counts optional matrix of raw experimental counts; if provided, will use to calculate baseline expression and guide-wise dispersion
#' @param effect_mode "fixed", "parametric", or "empirical": modes that determine how gene effects are simulated
#' @param gene_params optional list for parametric mode (e.g. list(mu = log(4), sd = 0.2))
#' @param guide_sd optional standard deviation for guide effect simulation
#' @param quantiles optional quantiles (e.g. c(0.05, 0.95)) to draw significant effects from in empirical mode
#' @return list with simulated counts, guide to gene mapping, sample design and true parameter values
#' @export
make_playdough <- function(n_genes,
                           guides_per_gene,
                           n_control = 2,
                           n_treatment = 2,
                           p_effects = 0.2,
                           p_positive = 0.1,
                           fold_change = 4,
                           n_ntc = 1000,
                           seed = 13,
                           counts = NULL,
                           effect_mode = c("fixed", "parametric", "empirical"),
                           gene_params = list(mu = log(fold_change), sd = 0.2),
                           guide_sd = 0.2,
                           quantiles = c(0.05, 0.95)) {
  set.seed(seed)
  effect_mode <- match.arg(effect_mode)
  n_samples  <- n_control + n_treatment
  n_effects  <- ceiling(p_effects * n_genes)
  n_guides <- n_genes * guides_per_gene + n_ntc

  # create guide to gene mapping
  guide_ids <- 1:n_guides
  gene_map <-
    c(rep(1:n_genes, each = guides_per_gene), rep(n_genes + 1, n_ntc))
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
    row_vars  <- matrixStats::rowVars(counts_mat)

    beta0_real <- log(row_means + 1)
    phi_real   <- pmax(1 / (row_vars / row_means ^ 2 - 1), 1e-3)

    # randomly sample pairs of moments (beta0, phi)
    moment_indices <-
      sample(nrow(counts_mat), n_guides, replace = TRUE)
    beta0_g <- beta0_real[moment_indices]
    phi_g <- phi_real[moment_indices]
  }
  else {
    beta0_g <- rnorm(n_guides, 5, 0.5)
    a <- runif(1, 0, 1)
    b <- runif(1, 0, 10)
    phi_g <- abs(a + (b / exp(beta0_g)) + rnorm(n_guides, 0, 0.05))
  }

  # simulate gene-level effects
  gene_effect <- numeric(n_genes)

  if (effect_mode == "fixed") {
    # fixed gene effect (one numeric value)
    signs <- sample(c(rep(-1, round(
      n_effects * (1 - p_positive)
    )),
    rep(1, n_effects - round(
      n_effects * (1 - p_positive)
    ))))
    gene_effect[1:n_effects] <- log(fold_change) * signs

  } else if (effect_mode == "parametric") {
    # parametric gene effect (distribution with given parameters)
    gene_mu <- gene_params$mu
    gene_sd <- gene_params$sd

    signs <- sample(c(rep(-1, round(
      n_effects * (1 - p_positive)
    )),
    rep(1, n_effects - round(
      n_effects * (1 - p_positive)
    ))))
    gene_effect[1:n_effects] <-
      abs(rnorm(n_effects, mean = gene_mu, sd = gene_sd)) * signs

  } else if (effect_mode == "empirical") {
    # empirical gene effect (distribution estimated from empirical data)
    lfc <-
      log(row_means + 1) - median(log(row_means) + 1) # empirical log-fold change, blind to sample conditions
    if (!is.null(quantiles)) {
      # take extreme effects based on defined quantiles
      q <- quantile(lfc, probs = quantiles, na.rm = TRUE)
      lfc <- lfc[(lfc <= q[1] | lfc >= q[2])]
    }

    gene_effect[1:n_effects] <-
      sample(lfc, n_effects, replace = TRUE)
  }

  # simulate guide-level effects
  guide_sd <- guide_sd
  ntc_sd <- 0.1
  beta1_g <- rnorm(length(gene_map[1:(n_genes * guides_per_gene)]),
                   mean = gene_effect[gene_map[1:(n_genes * guides_per_gene)]],
                   sd   = guide_sd)
  beta1_g <- c(beta1_g, rnorm(n_ntc, mean = 0, sd = ntc_sd))


  # calculate baseline mean counts
  mu_mat <- matrix(exp(beta0_g), nrow = n_guides, ncol = n_samples)

  # simulate sample-wise dispersion scalars
  gamma <- c(exp(rnorm(n_control, -0.5, 0.25)),
             exp(rnorm(n_treatment, 0.5, 0.25)))

  # simulate size factors
  size_factors <- runif(n_samples, 0.8, 1.2)

  # apply treatment effect
  treatment_indices <- (n_control + 1):n_samples
  targeting_indices <- 1:(n_genes * guides_per_gene)
  mu_mat[targeting_indices, treatment_indices] <-
    exp(beta0_g[targeting_indices] + beta1_g[targeting_indices])
  mu_mat <-
    mu_mat * matrix(size_factors,
                    nrow = n_guides,
                    ncol = n_samples,
                    byrow = TRUE)

  #  calculate dispersion
  size_mat <- 1 / (phi_g %o% gamma)

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
                              design = c(rep("control", n_control), rep("treatment", n_treatment)))

  controls <-
    data.frame(guides = guide_to_gene$sgRNA[guide_to_gene$gene == (n_genes + 1)],
               index = which(guide_to_gene$gene == (n_genes + 1)))

  playdough <- list(
    data = list(
      counts   = sim_counts,
      row_data = guide_to_gene,
      col_data = sample_design,
      controls = controls
    ),
    truth = list(
      significant_genes = which(gene_effect != 0),
      significant_guides = which(gene_map %in% which(gene_effect != 0)),
      beta0  = beta0_g,
      beta1  = beta1_g,
      mu     = gene_effect,
      phi    = phi_g,
      gamma  = gamma,
      sf = size_factors
    )
  )

  return(playdough)
}
