##### quality control plots #####

#' plot distribution of log2 (optional: normalized) counts per sample
#'
#' @param dough dough object with $data$counts
#' @param normalized logical indicating whether to normalize counts or not
#' @return p ggplot object
#' @export
plot_counts_violin <- function(dough, normalize = TRUE) {
  if (is.null(dough$data$counts))
    stop("no counts found")

  if (normalize) {
    norm_counts <- normalize_counts(dough)
    counts <- as.data.frame(norm_counts)
    title <- "distribution of log-normalized counts per sample"
    ylabel <- "log(normalized counts + 1)"
  }
  else {
    counts <- as.data.frame(dough$data$counts)
    title <- "distribution of log counts per sample"
    ylabel <- "log(counts + 1)"
  }

  if (all(grepl("^V[0-9]+$", colnames(counts)))) {
    colnames(counts) <- as.character(seq_len(ncol(counts)))
  }

  counts$guide <- rownames(counts)
  counts <- pivot_longer(counts,-guide, names_to = "sample", values_to = "counts")
  p <- ggplot(counts, aes(
      x = sample,
      y = log2(counts + 1),
      fill = sample
    )) +
    geom_violin(trim = FALSE, alpha = 0.5) +
    geom_boxplot(width = 0.1, outlier.size = 0.5) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = title, y = ylabel)

  return(p)
}

#' plot histogram of read counts per sample
#'
#' @param dough dough object with $data$counts
#' @param normalize logical indicating whether to normalize counts or not
#' @return p ggplot object
#' @export
plot_counts_density <- function(dough, normalize = TRUE) {
  if (is.null(dough$data$counts))
    stop("no counts found")

  if (normalize) {
    norm_counts <- normalize_counts(dough)
    counts <- as.data.frame(norm_counts)
    title <- "distribution of log-normalized counts per sample"
    xlabel <- "log2(normalized counts + 1)"
  } else {
    counts <- as.data.frame(dough$data$counts)
    title <- "distribution of log counts per sample"
    xlabel <- "log2(counts + 1)"
  }

  if (all(grepl("^V[0-9]+$", colnames(counts)))) {
    colnames(counts) <- as.character(seq_len(ncol(counts)))
  }

  counts$guide <- rownames(counts)
  counts <- pivot_longer(counts,-guide, names_to = "sample", values_to = "count")

  p <- ggplot(counts, aes(x = log2(count + 1), fill = sample)) +
    geom_density(alpha = 0.5) +
    labs(title = title, x = xlabel, y = "density") +
    theme_minimal()

  return(p)
}

#' plot sample correlation
#'
#' @param dough dough object with $counts
#' @return p ggplot object
#' @export
plot_sample_correlation <- function(dough) {
  if (is.null(dough$data$counts))
    stop("no counts found")

  norm_counts <- normalize_counts(dough)
  log_norm <- log2(norm_counts + 1)
  corr <- cor(norm_counts, method = "pearson")
  p <- pheatmap(corr, main = 'sample-wise log-normalized counts correlation')

  return(p)
}

#' plot moment estimated dispersion
#'
#' @param dough dough object with $data$counts
#' @return p ggplot object
#' @export
plot_dispersion <- function(dough) {
  if (is.null(dough$data$counts))
    stop("no counts found in dough$data$counts")

  counts <- as.matrix(dough$data$counts)

  # calculate mean and variance per guide
  guide_means <- rowMeans(counts)
  guide_vars  <- matrixStats::rowVars(counts)

  # method-of-moments NB dispersion
  # phi = (variance - mean) / mean^2
  phi_est <- pmax(((guide_vars - guide_means) / guide_means^2), 0)  # avoid negatives

  df <- data.frame(
    mean_count = guide_means,
    phi = phi_est,
    guide = dough$data$row_data$sgRNA,
    gene = dough$data$row_data$gene
  )

  # plot phi vs mean count
  p <- ggplot(df, aes(x = mean_count, y = phi)) +
    geom_point(alpha = 0.5) +
    scale_x_log10() +  # log scale for counts
    scale_y_log10() +  # log scale for dispersion
    theme_minimal() +
    xlab("mean counts per guide (log scale)") +
    ylab("estimated dispersion (log scale)") +
    ggtitle("method-of-moments estimate dispersion per guide")

  return(p)
}

#' plot MA of guides using empirical estimates (treatment vs control)
#'
#' @param dough dough object with $data$counts
#' @return p ggplot object
#' @export
plot_count_logfc <- function(dough) {
  if (is.null(dough$data$counts))
    stop("no counts found")

  sample_design <- dough$data$col_data$design
  guide_names <- dough$data$row_data$sgRNA

  if (is.null(sample_design)) {
    stop("col_data must have a 'sample_design' vector with values 'control' and 'treatment'")
  }

  # get sample indices
  control_samples <- which(sample_design == "control")
  treatment_samples <- which(sample_design == "treatment")

  if (length(control_samples) == 0 |
      length(treatment_samples) == 0) {
    stop("need both control and treatment samples in sample design")
  }

  # log-normalized counts with pseudocount
  norm_counts <- normalize_counts(dough)
  norm_counts <- log2(norm_counts + 1)

  # average across replicates
  control_mean <- rowMeans(norm_counts[, control_samples, drop = FALSE])
  treatment_mean <- rowMeans(norm_counts[, treatment_samples, drop = FALSE])

  # calculate log fold change
  logFC <- treatment_mean - control_mean
  mean_abundance <- (treatment_mean + control_mean) / 2

  ma <- data.frame(mean_abundance, logFC, guide = guide_names)

  p <- ggplot(ma, aes(x = mean_abundance, y = logFC)) +
    geom_point(alpha = 0.5) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    theme_minimal() +
    xlab("mean log2 normalized counts") +
    ylab("log2 fold change (treatment vs control)") +
    ggtitle("empirical guide log fold change vs abundance")

  return(p)
}

#' plot density of empirical logFC of guides, startified by non-targeting and targeting guides
#'
#' @param dough dough object with $data$counts
#' @return p ggplot object
#' @export
plot_logfc_density <- function(dough) {
  if (is.null(dough$data$counts))
    stop("no counts found")

  sample_design <- dough$data$col_data$design
  row_data <- dough$data$row_data

  # sample indices
  control_samples <- which(sample_design == "control")
  treatment_samples <- which(sample_design == "treatment")
  if (length(control_samples) == 0 || length(treatment_samples) == 0)
    stop("need both control and treatment samples in sample design")

  # log-normalized counts with pseudocount
  norm_counts <- normalize_counts(dough)
  norm_counts <- log2(norm_counts + 1)

  # calculate logFC
  control_mean <- rowMeans(norm_counts[, control_samples, drop = FALSE])
  treatment_mean <- rowMeans(norm_counts[, treatment_samples, drop = FALSE])
  logFC <- treatment_mean - control_mean

  # guide type
  guide_type <- ifelse(row_data$sgRNA %in% dough$data$controls$guide,
                       'non-targeting', 'targeting')

  df <- data.frame(logFC = logFC, type = guide_type)

  # density plot
  p <- ggplot(df, aes(x = logFC, fill = type)) +
    geom_density(alpha = 0.5) +
    theme_minimal() +
    xlab("log2 fold change (treatment vs control)") +
    ylab("density") +
    ggtitle("density of logFC for targeting vs non-targeting control guides")

  return(p)
}


#' plot histogram of number of guides per gene
#'
#' @param dough dough object with $data$row_data
#' @return ggplot object
#' @export
plot_guides_per_gene <- function(dough) {
  if (is.null(dough$data$row_data))
    stop("no guide-to-gene mapping found")

  guides <- as.data.frame(dough$data$row_data)

  # filter out non-targeting guides if ntc exists
  if (!is.null(dough$data$controls)) {
    guides <- guides[!guides$sgRNA %in% dough$data$controls$guide, , drop = FALSE]
  }

  guides_per_gene <- as.data.frame(table(guides$gene))
  colnames(guides_per_gene) <- c("gene", "n_guides")

  p <- ggplot(guides_per_gene, aes(x = n_guides)) +
    geom_histogram(binwidth = 1) +
    labs(title = "distribution of guides per gene",
         x = "# guides per gene",
         y = "# genes") +
    theme_minimal()

  return(p)
}

##### fit analysis plots #####

#' plot densities of targeting vs non-targeting guide log fold changes
#'
#' @param biscuit a fitted biscuit object with $results, $data
#' @return p ggplot object
#' @export
plot_guide_density <- function(biscuit) {
  if (is.null(biscuit$results$beta1))
    stop("no beta1 results found")

  # annotate sgRNAs as targeting or non-targeting
  beta1_summary <- biscuit$results$beta1 %>%
    mutate(control = if (!is.null(biscuit$data$controls)) {
      ifelse(sgRNA %in% biscuit$data$controls$guide,
             "non-targeting",
             "targeting")
    } else {
      "targeting"
    })

  # plot densities
  p <- ggplot(beta1_summary, aes(x = mean, fill = control)) +
    geom_density(alpha = 0.5) +
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
#' @param biscuit fitted biscuit object with $results$beta1 and $data$controls
#' @return ggplot object
#' @export
plot_guide_violin <- function(biscuit) {
  if (is.null(biscuit$results$beta1))
    stop("no beta1 results found")

  beta1_summary <- biscuit$results$beta1 %>%
    mutate(type = ifelse(
      sgRNA %in% biscuit$data$controls$guide,
      "non-targeting",
      "targeting"
    ))

  p <- ggplot(beta1_summary, aes(x = type, y = mean, fill = type)) +
    geom_violin(trim = FALSE, alpha = 0.5) +
    geom_boxplot(width = 0.1, outlier.shape = NA) +
    labs(x = NULL,
         y = "posterior mean guide effect (beta1)",
         fill = 'sgRNA type') +
    theme_minimal()

  return(p)
}


#' plot score versus rank plot of genes
#'
#' @param biscuit a biscuit object with $results
#' @param lfdr_threshold threshold to highlight significant genes
#' @param top_n number of top genes to label
#' @return ggplot object
#' @export
plot_gene_rank <- function(biscuit, lfdr_threshold = 0.05, top_n = 10) {
    if (is.null(biscuit$results$mu))
      stop("no mu results found")

    mu_summary <- biscuit$results$mu %>%
      arrange(desc(mean)) %>%
      mutate(
        rank = row_number(),
        category = case_when(
          lfdr < lfdr_threshold & mean > 0 ~ "positive",
          lfdr < lfdr_threshold & mean < 0 ~ "negative",
          TRUE ~ "not significant"
        )
      )

    # pick top_n significant genes with lowest lfdr
    top_genes <- mu_summary %>%
      filter(lfdr < lfdr_threshold) %>%
      arrange(desc(abs(mean))) %>%
      head(top_n)

    p <- ggplot(mu_summary, aes(x = rank, y = mean, color = category)) +
      geom_point(alpha = 0.5) +
      scale_color_manual(values = c(
        "not significant" = "grey",
        "positive" = "#F8766D",
        "negative" = "#00BFC4"
      )) +
      geom_text_repel(
        data = top_genes,
        aes(label = gene),
        nudge_y = 0.05 * max(mu_summary$mean),
        size = 3,
        segment.alpha = 0.5,
        max.overlaps = 20,
        show.legend = FALSE
      ) +
      labs(
        x = "gene, descending by posterior mean gene effect (mu)",
        y = "posterior mean gene effect (mu)",
        color = paste0("significance (lfdr < ", lfdr_threshold, ")")
      ) +
      theme_minimal()

    return(p)
  }

#'plot volcano plot for inferred gene-level effects
#'
#' @param biscuit a biscuit object with $results
#' @param lfdr_threshold threshold to highlight significant genes
#' @param top_n number of top genes to label
#' @return ggplot object
#' @export
plot_gene_volcano <- function(biscuit, lfdr_threshold = 0.05, top_n = 10) {
    if (is.null(biscuit$results$mu))
      stop("no mu results found")
    mu_summary <- biscuit$results$mu

    # flag significant genes
    mu_summary <- mu_summary %>%
      mutate(
        category = case_when(
          lfdr < lfdr_threshold & mean > 0 ~ "positive",
          lfdr < lfdr_threshold & mean < 0 ~ "negative",
          TRUE ~ "not significant"
        )
      )

    # select top_n genes by effect magnitude for labeling
    top_genes <- mu_summary %>%
      filter(lfdr < lfdr_threshold) %>%
      arrange(desc(abs(mean))) %>%
      head(top_n)

    # plot
    p <-
      ggplot(mu_summary, aes(
        x = mean,
        y = -log10(pmax(lfdr, 1e-5)),
        color = category
      )) +
      geom_point(alpha = 0.5) +
      coord_cartesian(clip = "off") +
      scale_color_manual(values = c(
        "not significant" = "grey",
        "positive" = "#F8766D",
        "negative" = "#00BFC4"
      )) +
      theme_minimal() +
      labs(
        x = "posterior mean gene effect (mu)",
        y = "-log10(lfdr)",
        color = paste0("significance (lfdr < ", lfdr_threshold, ")")
      ) +
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



#' plot MA plot of guides using inferred guide effects (beta1)
#'
#' @param biscuit biscuit object with $results
#' @param lfdr_threshold threshold to highlight significant guides
#' @param top_n number of top guides to label
#' @return ggplot object
#' @export
plot_guide_ma <- function(biscuit, lfdr_threshold = 0.05, top_n = 10) {
    if (is.null(biscuit$results$beta1))
      stop("no beta1 results found")

    norm_counts <- normalize_counts(biscuit)
    mean_abundance <- log2(rowMeans(norm_counts) + 1)

    # combine
    beta1_summary <- biscuit$results$beta1 %>%
      mutate(
        mean_abundance = mean_abundance[index],
        logFC = mean,
        category = case_when(
          lfdr < lfdr_threshold & mean > 0 ~ "positive",
          lfdr < lfdr_threshold & mean < 0 ~ "negative",
          TRUE ~ "not significant"
        )
      )

    # select top_n significant genes by effect magnitude for labeling
    top_guides <- beta1_summary %>%
      filter(lfdr < lfdr_threshold) %>%
      arrange(desc(abs(mean))) %>%
      head(top_n)

    # plot
    p <- ggplot(beta1_summary, aes(x = mean_abundance, y = logFC, color = category)) +
      geom_point(alpha = 0.5) +
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
      geom_hline(yintercept = 0,
                 linetype = "dashed",
                 color = "black") +
      labs(
        x = "mean(log2 normalized counts)",
        y = "posterior mean guide effect (beta1)",
        color = paste0("significance (lfdr < ", lfdr_threshold, ")"),
        title = "guide-level MA plot"
      ) +
      theme_minimal()

    return(p)
  }
