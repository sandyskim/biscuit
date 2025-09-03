##### helper functions #####

#' normalize counts
#'
#' @param dough dough (or biscuit) object with $data$counts
#' @return norm_counts matrix of normalized read counts, via median of ratios method, using DESeq2
get_norm_counts <- function(dough) {
  counts <- dough$data$counts
  sf <- estimateSizeFactorsForMatrix(counts)
  norm_counts <- t(t(counts) / sf)
  return(norm_counts)
}

##### quality control plots #####

#' plot distribution of log2 (optional: normalized) counts per sample
#'
#' @param dough dough object with $data$counts
#' @param normalized logical indicating whether to normalize counts or not
#' @return p ggplot object
plot_counts_boxplot <- function(dough, normalized = TRUE) {
  if (normalized) {
    norm_counts <- get_norm_counts(dough)
    df <- as.data.frame(norm_counts)
    title <- "distribution of log-normalized counts per sample"
    ylabel <- "log(normalized counts + 1)"
  }
  else {
    df <- as.data.frame(dough$data$counts)
    title <- "distribution of log counts per sample"
    ylabel <- "log(counts + 1)"
  }

  if (all(grepl("^V[0-9]+$", colnames(df)))) {
    colnames(df) <- as.character(seq_len(ncol(df)))
  }

  df$guide <- rownames(df)
  df <- pivot_longer(df, -guide, names_to = "sample", values_to = "counts")
  p <- ggplot(df, aes(x = sample, y = log2(counts + 1), fill = sample)) +
    geom_violin(trim = FALSE, alpha = 0.6) +
    geom_boxplot(width = 0.1, outlier.size = 0.5, alpha = 0.8) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = title, y = ylabel)

  return(p)
}

#' plot histogram of read counts per sample
#'
#' @param dough object with $data$counts
#' @param normalized logical indicating whether to normalize counts or not
#' @return p ggplot object
plot_counts_distribution <- function(dough, normalized = TRUE) {
  if (normalized) {
    norm_counts <- get_norm_counts(dough)
    df <- as.data.frame(norm_counts)
    title <- "distribution of log-normalized counts per sample"
    xlabel <- "log2(normalized counts + 1)"
  } else {
    df <- as.data.frame(dough$data$counts)
    title <- "distribution of log counts per sample"
    xlabel <- "log2(counts + 1)"
  }

  if (all(grepl("^V[0-9]+$", colnames(df)))) {
    colnames(df) <- as.character(seq_len(ncol(df)))
  }

  df$guide <- rownames(df)

  df <- pivot_longer(
    df,
    -guide,
    names_to = "sample",
    values_to = "count"
  )

  p <- ggplot(df, aes(x = log2(count + 1), fill = sample)) +
    geom_density(alpha = 0.5) +
    labs( title = title, x = xlabel, y = "frequency") +
    theme_minimal()

  return(p)
}

#' plot sample correlation
#'
#' @param dough dough object with $counts
#' @return p ggplot object
plot_sample_correlation <- function(dough) {
  norm_counts <- get_norm_counts(dough)
  log_norm <- log2(norm_counts + 1)
  corr <- cor(norm_counts, method = "pearson")
  p <- pheatmap(corr, main = 'sample-wise log-normalized counts correlation')

  return(p)
}

#' plot MA of guides using empirical estimates (treatment vs control)
#'
#' @param dough dough object with $counts
#' @return p ggplot object
plot_count_logfc <- function(dough) {
  sample_design <- dough$data$col_data$design
  guide_names <- dough$data$row_data$sgRNA

  if (is.null(sample_design)) {
    stop("col_data must have a 'sample_design' vector with values 'control' and 'treatment'")
  }

  # get sample indices
  control_samples <- which(sample_design == "control")
  treatment_samples <- which(sample_design == "treatment")

  if (length(control_samples) == 0 | length(treatment_samples) == 0) {
    stop("need both control and treatment samples in sample design")
  }

  # log-normalized counts with pseudocount
  norm_counts <- get_norm_counts(dough)
  norm_counts <- log2(norm_counts + 1)

  # average across replicates
  control_mean <- rowMeans(norm_counts[, control_samples, drop = FALSE])
  treatment_mean <- rowMeans(norm_counts[, treatment_samples, drop = FALSE])

  # calculate log fold change
  logFC <- treatment_mean - control_mean
  mean_abundance <- (treatment_mean + control_mean) / 2

  df <- data.frame(mean_abundance, logFC, guide = guide_names)

  p <- ggplot(df, aes(x = mean_abundance, y = logFC)) +
    geom_point(alpha = 0.4) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    theme_minimal() +
    xlab("mean log2 normalized counts") +
    ylab("log2 fold change (treatment vs control)") +
    ggtitle("empirical guide log fold change vs abundance")

  return(p)
}

#' plot histogram of number of guides per gene
#'
#' @param dough dough object with $data$row_data
#' @return ggplot object
plot_guides_per_gene <- function(dough) {
  df <- as.data.frame(dough$data$row_data)

  # filter out non-targeting guides only if ntc exists
  if (!is.null(dough$data$ntc)) {
    df <- df %>% filter(!(sgRNA %in% dough$data$ntc$guide))
  }

  df <- df %>% count(gene, name = "n_guides")

  p <- ggplot(df, aes(x = n_guides)) +
    geom_histogram(binwidth = 1, fill = "steelblue", color = "white") +
    labs(title = "Distribution of sgRNAs per gene",
         x = "# sgRNAs per gene", y = "# genes") +
    theme_minimal()

  return(p)
}

##### fit analysis plots #####

#' plot densities of targeting vs non-targeting guide log fold changes
#'
#' @param biscuit a fitted biscuit object with $results$beta1 and $data$ntc
#' @return ggplot object
plot_guide_density <- function(biscuit) {
  # annotate sgRNAs as targeting or non-targeting
  beta1_summary <- biscuit$results$beta1 %>%
    mutate(control = if(!is.null(biscuit$data$ntc)) {
      ifelse(sgRNA %in% biscuit$data$ntc$guide, "non-targeting", "targeting")
    } else {
      "targeting"
    })

  # plot densities
  p <- ggplot(beta1_summary, aes(x = mean, fill = control)) +
    geom_density(alpha = 0.4) +
    geom_vline(xintercept = 0, linetype = "dotted") +
    labs(
      title = "distribution of log fold changes",
      x = "posterior mean log fold change",
      y = "density",
      fill = "sgRNA type"
    ) +
    theme_minimal()

  return(p)
}


#' plot violin plot of biscuit targeting vs non-targeting control guide log fold changes
#'
#' @param biscuit fitted biscuit object with $results$beta1 and $data$ntc
#' @return ggplot object
plot_guide_violin <- function(biscuit) {
  if (is.null(biscuit$results$beta1)) stop("no beta1 results found in biscuit$results")

  beta1_summary <- biscuit$results$beta1 %>%
    mutate(
      type = ifelse(
        sgRNA %in% biscuit$data$ntc$guide,
        "non-targeting",
        "targeting"
      )
    )

  ggplot(beta1_summary, aes(x = type, y = mean, fill = type)) +
    geom_violin(trim = FALSE, alpha = 0.6) +
    geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.8) +
    labs(
      x = NULL,
      y = "posterior mean guide effect (beta1)",
      fill = 'sgRNA type'
    ) +
    theme_minimal()
}


#' plot score versus rank plot of genes
#'
#' @param biscuit a biscuit object with $results
#' @param lfsr_threshold threshold to highlight significant genes
#' @param top_n number of top genes to label
#' @return ggplot object
plot_gene_rank <- function(biscuit, lfsr_threshold = 0.05, top_n = 10) {
  if (is.null(biscuit$results$mu)) stop("no mu results found in biscuit$results")

  mu_summary <- biscuit$results$mu %>%
    arrange(desc(mean)) %>%
    mutate(
      rank = row_number(),
      category = case_when(
        lfsr < lfsr_threshold & mean > 0 ~ "positive",
        lfsr < lfsr_threshold & mean < 0 ~ "negative",
        TRUE ~ "not significant"
      )
    )

  # pick top_n genes with lowest lfsr
  top_genes <- mu_summary %>%
    filter(lfsr < lfsr_threshold) %>%
    arrange(desc(abs(mean))) %>%
    head(top_n)

  ggplot(mu_summary, aes(x = rank, y = mean, color = category)) +
    geom_point(alpha = 0.7) +
    scale_color_manual(
      values = c(
        "not significant" = "grey",
        "positive" = "#F8766D",
        "negative" = "#00BFC4"
      )
    ) +
    geom_text_repel(
      data = top_genes, aes(label = gene),
      nudge_y = 0.05 * max(mu_summary$mean),
      size = 3, segment.alpha = 0.5, max.overlaps = 20,
      show.legend = FALSE
    ) +
    labs(
      x = "gene rank",
      y = "posterior mean gene effect (mu)",
      color = paste0("significance (lfsr < ", lfsr_threshold, ")")
    ) +
    theme_minimal()
}

#'plot volcano plot for inferred gene-level effects
#'
#' @param biscuit a biscuit object with $results
#' @param lfsr_threshold threshold to highlight significant genes
#' @param top_n number of top genes to label
#' @return ggplot object
plot_gene_volcano <- function(biscuit, lfsr_threshold = 0.05, top_n = 10) {
  mu_summary <- biscuit$results$mu

  # flag significant genes
  mu_summary <- mu_summary %>%
    mutate(category = case_when(
      lfsr < lfsr_threshold & mean > 0 ~ "positive",
      lfsr < lfsr_threshold & mean < 0 ~ "negative",
      TRUE ~ "not significant"
    ))

  # select top_n genes by effect magnitude for labeling
  top_genes <- mu_summary %>%
    filter(lfsr < lfsr_threshold) %>%
    arrange(desc(abs(mean))) %>%
    head(top_n)

  # plot
  p <- ggplot(mu_summary, aes(x = mean, y = -log10(pmax(lfsr, 1e-5)), color = category)) +
    geom_point(alpha = 0.8) +
    coord_cartesian(clip = "off") +
    scale_color_manual(values = c(
      "not significant" = "grey",
      "positive" = "#F8766D",
      "negative" = "#00BFC4"
    )) +
    theme_minimal() +
    labs(x = "posterior mean gene effect (mu)", y = "-log10(lfsr)",
         color = paste0("significance (lfsr < ", lfsr_threshold, ")")) +
    geom_text_repel(
      data = top_genes,
      aes(label = gene),
      size = 3,
      max.overlaps = 20,
      show.legend = FALSE
    ) +
    theme(legend.position = "right")

  return(p)
}

#' plot stacked plot of posterior densities of a given gene and the guides that target it
#'
#' @param biscuit a biscuit object with $fit
#' @param gene_name name of gene
#' @return ggplot object

plot_mu_beta1_density <- function(biscuit, gene_name) {
  if (is.null(biscuit$results$mu)) stop("no mu results found in biscuit$results")
  if (is.null(biscuit$results$beta1)) stop("no beta1 results found in biscuit$results")

  # find gene index in row_data
  unique_genes <- unique(biscuit$data$row_data[,2])
  gene_index <- which(unique_genes == gene_name)
  if (length(gene_index) == 0) stop("gene not found in biscuit$data$row_data")

  # extract posterior draws for the gene effect
  mu_col <- paste0("mu[", gene_index, "]")
  mu_draws <- as.data.frame(biscuit$fit$posterior)[, mu_col, drop = FALSE] %>%
    as_tibble() %>%
    pivot_longer(cols = everything(), names_to = "parameter", values_to = "value") %>%
    mutate(label = gene_name,
           type = "gene")

  # extract posterior draws for all guides targeting this gene
  guides_idx <- which(biscuit$data$row_data[,2] == gene_name)
  guide_names <- biscuit$data$row_data[guides_idx, 1]

  guide_draws <- lapply(seq_along(guides_idx), function(i) {
    idx <- guides_idx[i]
    col <- paste0("beta1[", idx, "]")
    as.data.frame(biscuit$fit$posterior)[, col, drop = FALSE] %>%
      as_tibble() %>%
      pivot_longer(cols = everything(), names_to = "parameter", values_to = "value") %>%
      mutate(label = guide_names[i],
             type = "guide")
  }) %>% bind_rows()

  # combine posterior draws of mu and corresponding beta1
  plot_df <- bind_rows(mu_draws, guide_draws)

  # reorder so mu is first
  plot_df$label <- factor(plot_df$label, levels = rev(c(gene_name, guide_names)))

  # plot stacked densities
  ggplot(plot_df, aes(x = value, y = label, fill = type)) +
    ggridges::geom_density_ridges(scale = 1.5, alpha = 0.7) +
    labs(
      x = "posterior draws",
      y = NULL,
      title = paste0("posterior densities for gene ", gene_name, " and its guides"),
      fill = NULL
    ) +
    theme_minimal() +
    theme(legend.position = "none")
}

#' plot stacked density plot of phi, and corresponding phi times gamma, visualizing variation in dispersion across samples
#'
#' @param biscuit a biscuit object with $results and $fit
#' @return ggplot object

plot_phi_gamma_density <- function(biscuit) {
  if (is.null(biscuit$results$phi)) stop("no phi results found in biscuit$results")
  if (is.null(biscuit$results$gamma)) stop("no gamma results found in biscuit$results")

  # phi and gamma means
  phi_means <- biscuit$results$phi$mean

  # gamma means
  gamma_means <- biscuit$results$gamma$mean
  n_samples <- length(gamma_means)

  # phi
  phi <- data.frame(
    value = phi_means,
    label = "phi",
    type = "phi"
  )

  # phi * gamma[sample]
  phi_gamma_list <- lapply(seq_len(n_samples), function(i) {
    data.frame(
      value = phi_means * gamma_means[i],
      label = paste0("phi * gamma[", i, "]"),
      type = "phi * gamma"
    )
  })

  phis <- bind_rows(phi, phi_gamma_list)

  # order factor so phi is on top
  phis$label <- factor(phis$label, levels = rev(unique(phis$label)))

  ggplot(phis, aes(x = log(value), y = label, fill = type)) +
    ggridges::geom_density_ridges(scale = 1.5, alpha = 0.7) +
    labs(
      x = "posterior means",
      y = NULL,
      title = "guide-level phi and phi * gamma across samples",
      fill = NULL
    ) +
    theme_minimal() +
    theme(legend.position = "none")
}

#' plot density plot of posterior distribution of a given given and posterior means of guides that target it
#'
#' @param biscuit a biscuit object with $results and $fit
#' @param gene_name name of gene to look at
#' @param lfsr_threshold threshold to highlight significant gene
#' @return ggplot object
plot_mu_beta1 <- function(biscuit, gene_name, lfsr_threshold=0.05) {
  if (is.null(biscuit$results$mu)) stop("no mu results found in biscuit$results")
  if (is.null(biscuit$results$beta1)) stop("no beta1 results found in biscuit$results")

  # guides targeting the gene
  guides_idx <- which(biscuit$data$row_data$gene == gene_name)
  guide_names <- biscuit$data$row_data$guide[guides_idx]
  if (length(guides_idx) == 0) stop("gene not found in biscuit$data$row_data")

  # posterior draws for mu
  unique_genes <- unique(biscuit$data$row_data[,2])
  gene_index <- which(unique_genes == gene_name)
  mu_col <- paste0("mu[", gene_index, "]")
  mu_draws <- as.data.frame(biscuit$fit$posterior)[, mu_col, drop = FALSE] %>%
    as_tibble() %>%
    pivot_longer(cols = everything(), names_to = "parameter", values_to = "value")

  mu_mean <- mean(mu_draws$value)

  mu_draws <- mu_draws %>%
    mutate(index = as.integer(stringr::str_extract(parameter, "\\d+"))) %>%
    left_join(
      mu_summary %>% select(index, lfsr),
      by = "index"
    ) %>%
    mutate(
      category = case_when(
        value > 0 & lfsr < lfsr_threshold ~ "positive",
        value < 0 & lfsr < lfsr_threshold ~ "negative",
        TRUE ~ "not significant"
      )
    )

  beta1_summary <- biscuit$results$beta1 %>%
    filter(gene == gene_name) %>%
    arrange(match(sgRNA, guide_names)) %>%
    pull(mean)

  p1 <- ggplot(mu_draws, aes(x = value)) +
    geom_density(aes(fill = category), alpha = 0.5) +
    scale_fill_manual(values = c(
      "not significant" = "grey",
      "positive" = "#F8766D",
      "negative" = "#00BFC4"
    )) +
    geom_vline(xintercept = mu_mean, linetype = "dashed", size = 0.5) +
    labs(title = paste0("posterior density for ", gene_name, " and its guides' posterior means"),
         x = NULL,
         y = "density",
         fill = paste0("significance (lfsr < ", lfsr_threshold, ")")) +
    theme_minimal() +
    theme(axis.text = element_blank(),
          axis.ticks = element_blank())

  xmin <- min(c(mu_draws$value, beta1_summary))
  xmax <- max(c(mu_draws$value, beta1_summary))
  p2 <- ggplot() +
    geom_rect(aes(xmin = xmin, xmax = xmax, ymin = 0, ymax = 1), fill = "grey", alpha = 0.5) +
    geom_vline(xintercept = mu_mean, size = 2) +
    geom_vline(xintercept = beta1_summary, size = 0.5) +
    scale_y_continuous(breaks = 0.5, labels = gene_name, limits = c(0,1)) +
    coord_cartesian(xlim = c(xmin, xmax)) +
    labs(x = "effect size", y = NULL) +
    theme_minimal() +
    theme(
      axis.ticks = element_blank(),
      axis.title.y = element_blank(),
      panel.grid = element_blank()
    )

  p1 / p2 + plot_layout(heights = c(1, 0.1))
}

#' plot MA plot of guides using inferred guide effects (beta1)
#'
#' @param biscuit biscuit object with $results
#' @param lfsr_threshold threshold to highlight significant guides
#' @param top_n number of top guides to label
#' @return ggplot object
plot_beta1_ma <- function(biscuit, lfsr_threshold = 0.05, top_n = 10) {
  if (is.null(biscuit$results$beta1)) stop("No beta1 results found in biscuit$results")

  # get beta1 summary
  beta1_summary <- biscuit$results$beta1

  # normalized counts per guide
  norm_counts <- get_norm_counts(biscuit)

  # compute average abundance per guide
  avg_abundance <- log2(rowMeans(norm_counts) + 1)

  # combine
  beta1_summary <- beta1_summary %>%
    mutate(
      A = avg_abundance[index],           # average log2 counts
      M = mean,
      category = case_when(
      lfsr < lfsr_threshold & mean > 0 ~ "positive",
      lfsr < lfsr_threshold & mean < 0 ~ "negative",
      TRUE ~ "not significant"
  ))

  # select top_n genes by effect magnitude for labeling
  top_guides <- beta1_summary %>%
    filter(lfsr < lfsr_threshold) %>%
    arrange(desc(abs(mean))) %>%
    head(top_n)

  # plot
  ggplot(beta1_summary, aes(x = A, y = M, color = category)) +
    geom_point(alpha = 0.7, size = 1.5) +
    scale_color_manual(values = c(
      "not significant" = "grey",
      "positive" = "#F8766D",
      "negative" = "#00BFC4"
    )) +
    geom_text_repel(
      data = top_guides,
      aes(label = sgRNA),
      size = 3,
      max.overlaps = 20,
      show.legend = FALSE
    ) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    labs(
      x = "mean(log2 normalized counts)",
      y = "posterior mean guide effect (beta1)",
      color = paste0("significance (lfsr < ", lfsr_threshold, ")"),
      title = "guide-level MA plot"
    ) +
    theme_minimal()
}
