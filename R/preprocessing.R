#' make dough object from counts, guide_to_gene mapping, and design of samples (optional: list of names of non-targeting controls)
#'
#' @param counts n_guides x n_samples matrix holding the raw guide/sgRNA counts in each sample
#' @param guide_to_gene n_guides x 2 matrix holding the name of the guide/sgRNA in the first column and the name of the gene it's targeting in the second column.
#' @param sample_design n_samples x 2 holding the name of the samples (which match the column names of the counts) in the first column and the experimental design in the second column. the first sample must be a control sample.
#' @param controls (optional) vector containing the names of the non-targeting guides, all of which must be in the first column of guide_to_gene
#' @return dough object, with data stored in $data
#' @export
make_dough <- function(counts, guide_to_gene, sample_design, controls = NULL) {
  if (length(unique(sample_design)) != 2) {
    stop("sample design must have exactly two conditions.")
  }
  if (nrow(counts) != nrow(guide_to_gene)) {
    stop("dimension mismatch! counts has ", nrow(counts), " rows and guide_to_gene mapping has ", nrow(guide_to_gene), " rows.\n",
         "length of guide-to-gene mapping must equal the number of guides (rows of counts).")
  }
  if(ncol(counts) != nrow(sample_design)){
    stop("dimension mismatch! counts has ", ncol(counts), " columns and sample_design has ", nrow(sample_design), " columns.\n",
         "length of sample to condition mapping must equal the number of samples (columns of counts).")
  }
  if(!all.equal(colnames(counts), sample_design[,1])) {
    stop("sample names in the sample-to-condition mapping do not align with the column names of the counts matrix.\n",
         "each column of the counts matrix must correspond to a sample listed in the mapping, in the same order")
  }
  if(!all.equal(rownames(counts), guide_to_gene[,1])) {
    stop("guide names in the guide-to-gene mapping do not align with the row names of the counts matrix.\n",
        "each row of the counts matrix must correspond to a guide listed in the mapping, in the same order")
  }

  # reorder targeting guides first, non-targeting controls last
  if (!is.null(controls)) {
    targeting_idx <- which(!guide_to_gene[,1] %in% controls)
    ntc_idx       <- which(guide_to_gene[,1] %in% controls)

    new_order <- c(targeting_idx, ntc_idx)
    counts    <- counts[new_order, , drop = FALSE]
    guide_to_gene  <- guide_to_gene[new_order, , drop = FALSE]

    # record ntc guides and their row indices in reordered counts
    controls <- data.frame(
      guide = guide_to_gene[ntc_idx, 1],
      index = match(ntc_idx, new_order),
      stringsAsFactors = FALSE
    )
  }

  guide_to_gene <- as.data.frame(guide_to_gene)
  colnames(guide_to_gene) <- c('sgRNA', 'gene')

  sample_design <- data.frame(
    sample = sample_design[,1],
    design = as.factor(sample_design[,2]),
    stringsAsFactors = FALSE
  )

  dough <- list(
    data = list(
      counts   = as.matrix(counts),
      row_data = as.data.frame(guide_to_gene),
      col_data = sample_design,
      controls = controls
    )
  )

  class(dough) <- "dough"
  return(dough)
}

#' filter targeting guide counts
#'
#' @param dough dough object with $data: $data$counts, $data$row_data, $data$col_data (optional: $data$controls)
#' @param min_per_sample minimum total counts per guide across a min_prop samples to keep
#' @param min_prop minimum proportion of samples for min_per_sample to meet min_per_sample, number of samples is rounded to nearest integer
#' @param min_guides_per_gene minimum number of guides per gene to keep the gene
#' @param verbose logical to print out dimensions before and after filtering
#' @return dough object with filtered data
#' @export
trim_dough <- function(dough, min_per_sample = 1, min_prop = 0.2, min_guides_per_gene = 2, verbose = TRUE) {

  counts <- dough$data$counts
  row_data <- dough$data$row_data
  col_data <- dough$data$col_data
  controls <- dough$data$controls

  # indicate which guides are non-targeting controls
  is_ntc <- row_data$sgRNA %in% controls$guide

  # filtering all zero counts
  keep_zero <- rowSums(counts) > 0

  # filter targeting guides by minimum count across a proportion of samples, rounded to the nearest integer
  keep_counts <- rowSums(counts >= min_per_sample) >= round(min_prop * ncol(counts))

  # filter genes by minimum guide per gene
  gene_table <- table(row_data$gene)
  keep_genes <- row_data$gene %in% names(gene_table[gene_table >= min_guides_per_gene])

  # keep guides that satisfy counts and gene rules, or are NT controls
  keep <- (keep_counts & keep_zero & keep_genes) | is_ntc

  counts_filtered <- counts[keep, , drop = FALSE]
  row_data_filtered <- row_data[keep, , drop = FALSE]

  # reconstruct count matrix so targeting guides come first and then non-targeting guides
  if(!is.null(controls)) {
    is_ntc_filtered <- row_data_filtered$sgRNA %in% controls$guide
    counts_filtered <- rbind(counts_filtered[!is_ntc_filtered, , drop = FALSE],
                             counts_filtered[is_ntc_filtered, , drop = FALSE])
    row_data_filtered <- rbind(row_data_filtered[!is_ntc_filtered, , drop = FALSE],
                          row_data_filtered[is_ntc_filtered, , drop = FALSE])
    controls <- controls[controls$guide %in% row_data_filtered$sgRNA, , drop = FALSE]
    controls$index <- match(controls$guide, row_data_filtered$sgRNA)
  }

  # save filtered data to dough
  dough$data <- list(
    counts = counts_filtered,
    row_data = row_data_filtered,
    col_data = col_data,
    controls = controls
  )

  # print before filtering and after filtering
  if (verbose) {
    message("guides before: ", nrow(row_data),
            "; after: ", nrow(row_data_filtered))
    message("genes before: ", length(unique(row_data$gene)),
            "; after: ", length(unique(row_data_filtered$gene)))
    message("non-targeting controls preserved: ",
            ifelse(!is.null(nrow(controls)), nrow(controls), 'no non-targeting guides were provided'))
  }

  return(dough)
}


