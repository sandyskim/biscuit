#' simulate count data using moments estimated from real data
#'
#' @param counts matrix of raw experimental counts
#' @param n_genes number of genes to simulate
#' @param guides_per_gene number of guides per gene to simulate
#' @param n_control number of control samples
#' @param n_treatment number of treatment samples
#' @param p_effects proportion of effects
#' @param p_positive proportion of effects that are positive
#' @param fold_change absolute fold change of significant effects
#' @param n_ntc number of non-targeting control guides
#' @param seed numeric seed for reproducibility
#' @return list including simulated counts, guide to gene mapping, sample design and true parameter values
#' @export
make_playdough <- function(counts, n_genes, guides_per_gene,
                                                  n_control = 2, n_treatment = 2,
                                                  p_effects = 0.1, p_positive = 0.1,
                                                  fold_change = 4, n_ntc = 1000, seed = 13) {
  set.seed(seed)

  counts_mat <- as.matrix(counts)
  n_samples <- n_control + n_treatment
  n_effects <- ceiling(p_effects * n_genes)

  # method of moments from real data
  row_means <- rowMeans(counts_mat + 1)
  row_vars <- rowVars(counts_mat) # can use matrixStats::rowVars for speed
  beta0_real <- log(row_means)
  phi_real <- pmax(1 / (row_vars / row_means^2 - 1), 1e-3)

  # randomly sample moments (beta0, phi)
  moment_indices <- sample(nrow(counts_mat), n_guides, replace = TRUE)
  beta0_g <- beta0_real[moment_indices]
  phi_g <- phi_real[moment_indices]

  # create guide to gene mapping
  n_guides <- n_genes * guides_per_gene + n_ntc
  guide_ids <- 1:n_guides
  gene_map <- c(rep(1:n_genes, each = guides_per_gene), rep(n_genes + 1, n_ntc))
  guide_to_gene <- data.frame(sgRNA = guide_ids, gene = gene_map)

  # simulate gene-level effects
  signs <- sample(c(rep(-1, round(n_effects * (1 - p_positive))),
                    rep(1, n_effects - round(n_effects * (1 - p_positive)))))
  gene_effect <- numeric(n_genes)
  gene_effect[1:n_effects] <- log(fold_change) * signs

  # simulate guide-level effects
  beta1_g <- gene_effect[gene_map[1:(n_genes * guides_per_gene)]]
  beta1_g <- c(beta1_g, rep(0, n_ntc))  # non-targeting guides

  # Simulate baseline mean counts
  mu_mat <- matrix(exp(beta0_g), nrow = n_guides, ncol = n_samples)

  # simulate sample-wise dispersion scalars
  gamma <- exp(c(rnorm(n_control, 0.75, 0.25), rnorm(n_treatment, 1.25, 0.25)))

  # simulate size factors
  size_factors <- runif(n_samples, 0.8, 1.2)

  # apply treatment effect
  treatment_idx <- (n_control + 1):n_samples
  targeting_idx <- 1:(n_genes * guides_per_gene)
  mu_mat[targeting_idx, treatment_idx] <- exp(beta0_g[targeting_idx] + beta1_g[targeting_idx])
  mu_mat <- mu_mat * matrix(size_factors, nrow = n_guides, ncol = n_samples, byrow = TRUE)

  size_vec <- 1 / (phi_g * gamma)

  # negative binomial simulation
  sim_counts <- matrix(rnbinom(n_guides * n_samples,
                               mu = as.vector(mu_mat),
                               size = rep(size_vec, each = n_guides)),
                       nrow = n_guides, ncol = n_samples)

  rownames(sim_counts) <- guide_ids
  colnames(sim_counts) <- c(paste0("control", 1:n_control),
                            paste0("treatment", 1:n_treatment))

  sample_design <- data.frame(sample = colnames(sim_counts),
                              design = c(rep("control", n_control), rep("treatment", n_treatment)))

  ntcs <- data.frame(guide = paste0("ntc", 1:n_ntc), index = (n_genes + 1))

  playdough <- list(
    counts = sim_counts,
    row_data = guide_to_gene,
    col_data = sample_design,
    gene_effect = gene_effect,
    ntc = ntcs,
    truth = list(
      significant_genes = 1:n_effects,
      beta0 = beta0_g,
      beta1 = beta1_g,
      mu = gene_effect,
      phi = phi_g,
      gamma = gamma
    )
  )

  return(playdough)
}
