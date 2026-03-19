#' Create a dough object from counts, guide-to-gene mapping, and sample design
#'
#' @param counts n_guides x n_samples matrix of raw guide/sgRNA counts
#' @param guide_to_gene n_guides x 2 data frame with guide names in column 1 and target gene names in column 2
#' @param sample_design n_samples x 2 data frame with sample names in column 1 and condition in column 2;
#'   column names must match \code{colnames(counts)}; first condition listed is treated as control
#' @param controls optional character vector of non-targeting guide names; must appear in column 1 of \code{guide_to_gene}
#' @return dough object with data stored in \code{$data}
#' @export
make_dough <- function(counts, guide_to_gene, sample_design, controls = NULL) {
  # input validation
  if (length(unique(sample_design[, 2])) < 2)
    stop("sample_design must contain at least two conditions.")

  if (nrow(counts) != nrow(guide_to_gene))
    stop(sprintf(
      "Dimension mismatch: counts has %d rows but guide_to_gene has %d rows.",
      nrow(counts), nrow(guide_to_gene)
    ))

  if (ncol(counts) != nrow(sample_design))
    stop(sprintf(
      "Dimension mismatch: counts has %d columns but sample_design has %d rows.",
      ncol(counts), nrow(sample_design)
    ))

  if (!isTRUE(all.equal(colnames(counts), as.character(sample_design[, 1]))))
    stop("Column names of counts must match sample names in sample_design, in the same order.")

  if (!isTRUE(all.equal(rownames(counts), as.character(guide_to_gene[, 1]))))
    stop("Row names of counts must match guide names in guide_to_gene, in the same order.")

  # reorder: targeting guides first, NTCs last
  if (!is.null(controls)) {
    targeting_idx <- which(!guide_to_gene[, 1] %in% controls)
    ntc_idx <- which(guide_to_gene[, 1] %in% controls)
    new_order <- c(targeting_idx, ntc_idx)

    counts <- counts[new_order, , drop = FALSE]
    guide_to_gene <- guide_to_gene[new_order, , drop = FALSE]

    ntc_rows <- (length(targeting_idx) + 1):nrow(guide_to_gene)
    controls <- data.frame(
      guide = guide_to_gene[ntc_rows, 1],
      index = ntc_rows,
      stringsAsFactors = FALSE
    )
  }

  guide_to_gene <- as.data.frame(guide_to_gene)
  colnames(guide_to_gene) <- c("sgRNA", "gene")

  condition <- as.integer(as.factor(sample_design[, 2]))
  sample_design <- data.frame(
    sample = sample_design[, 1],
    design = as.factor(ifelse(condition == 1, "control", "treatment")),
    stringsAsFactors = FALSE
  )

  list(data = list(
      counts   = as.matrix(counts),
      row_data = guide_to_gene,
      col_data = sample_design,
      controls = controls
  ))
}


#' Filter guides and genes from a dough object by minimum count thresholds
#'
#' @param dough dough object with \code{$data}
#' @param min_per_sample minimum counts per guide in at least \code{min_prop} of samples
#' @param min_prop minimum proportion of samples that must meet \code{min_per_sample}
#' @param min_guides_per_gene minimum number of guides required to retain a gene
#' @param verbose logical; whether to print guide and gene counts before and after filtering
#' @return dough object with filtered \code{$data}
#' @export
trim_dough <- function(dough,
                       min_per_sample = 10,
                       min_prop = 0.2,
                       min_guides_per_gene = 2,
                       verbose = TRUE) {
  counts <- dough$data$counts
  row_data <- dough$data$row_data
  col_data <- dough$data$col_data
  controls <- dough$data$controls

  is_ntc <- row_data$sgRNA %in% controls$guide

  keep_nonzero <- rowSums(counts) > 0
  keep_counts <- rowSums(counts >= min_per_sample) >= round(min_prop * ncol(counts))
  keep_genes <- row_data$gene %in% names(which(table(row_data$gene) >= min_guides_per_gene))
  keep <- (keep_nonzero & keep_counts & keep_genes) | is_ntc

  counts_f <- counts[keep, , drop = FALSE]
  row_data_f <- row_data[keep, , drop = FALSE]

  # ensure targeting guides precede NTCs
  if (!is.null(controls)) {
    is_ntc_f <- row_data_f$sgRNA %in% controls$guide
    counts_f <- rbind(counts_f[!is_ntc_f, , drop = FALSE],
                        counts_f[is_ntc_f,  , drop = FALSE])
    row_data_f <- rbind(row_data_f[!is_ntc_f, , drop = FALSE],
                        row_data_f[is_ntc_f,  , drop = FALSE])

    controls <- controls[controls$guide %in% row_data_f$sgRNA, , drop = FALSE]
    controls$index <- match(controls$guide, row_data_f$sgRNA)
  }

  dough$data <- list(
    counts = counts_f,
    row_data = row_data_f,
    col_data = col_data,
    controls = controls
  )

  if (verbose) {
    message("Filtering finished!")
    message(sprintf("Guides: %d -> %d", nrow(row_data), nrow(row_data_f)))
    message(sprintf("Genes:  %d -> %d", length(unique(row_data$gene)), length(unique(row_data_f$gene))))
    message("Non-targeting controls preserved: ",
            if (!is.null(controls)) nrow(controls) else "none provided")
  }

  dough
}
