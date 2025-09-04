#' make dough object from counts, guide_to_gene mapping, and design of samples (optional: list of names of non-targeting controls)
#'
#' @param counts n_guides x n_samples matrix holding the raw guide/sgRNA counts in each sample
#' @param guide_to_gene n_guides x 2 matrix holding the name of the guide/sgRNA in the first column and the name of the gene it's targeting in the second column
#' @param sample_design vector of length n_samples indicating experimental condition, where the first entry must be a control
#' @param ntcs (optional) vector containing the names of the non-targeting guides, all of which must be in the first column of guide_to_gene
#' @return dough object, with data stored in $data
#' @export
make_dough <- function(counts, guide_to_gene, sample_design, ntcs = NULL) {
  if (length(unique(sample_design)) != 2) {
    stop("sample design must have exactly two conditions.")
  }
  if (nrow(counts) != nrow(guide_to_gene)) {
    stop("dimension mismatch! counts has ", nrow(counts), " rows and guide_to_gene mapping has ", nrow(guide_to_gene), " rows.\n",
         "length of guide to gene mapping must equal the number of guides (rows of counts).")
  }
  if(ncol(counts) != length(sample_design)){
    stop("dimension mismatch! counts has ", ncol(counts), " columns and sample_design has ", ncol(sample_design), " columns.\n",
         "length of sample design vector must equal the number of samples (columns of counts).")
  }

  # reorder targeting guides first, non-targeting controls last
  if (!is.null(ntcs)) {
    targeting_idx <- which(!guide_to_gene[,1] %in% ntcs)
    ntc_idx       <- which(guide_to_gene[,1] %in% ntcs)

    new_order <- c(targeting_idx, ntc_idx)
    counts    <- counts[new_order, , drop = FALSE]
    guide_to_gene  <- guide_to_gene[new_order, , drop = FALSE]

    # record ntc guides and their row indices in reordered counts
    ntcs <- data.frame(
      guide = guide_to_gene[ntc_idx, 1],
      index = match(ntc_idx, new_order),
      stringsAsFactors = FALSE
    )
  }

  guide_to_gene <- as.data.frame(guide_to_gene)
  colnames(guide_to_gene) <- c('sgRNA', 'gene')

  sample_names <- colnames(counts)
  if (is.null(sample_names)) {
    sample_names <- paste0("sample", seq_len(ncol(counts)))
    colnames(counts) <- sample_names
  }
  sample_design <- data.frame(
    sample = sample_names,
    design = sample_design,
    stringsAsFactors = FALSE
  )

  dough <- list(
    data = list(
      counts   = as.matrix(counts),
      row_data = as.data.frame(guide_to_gene),
      col_data = sample_design,
      ntc = ntcs
    )
  )

  class(dough) <- "dough"
  return(dough)
}

#' filter targeting guide counts
#'
#' @param dough dough object with $data: $data$counts, $data$row_data, $data$col_data (optional: $data$ntc)
#' @param min_counts minimum total counts per guide across samples to keep
#' @param min_guides_per_gene minimum number of guides per gene to keep the gene
#' @param verbose logical to print out dimensions before and after filtering
#' @return dough object with filtered data
#' @export
trim_dough <- function(dough, min_counts = 30, min_guides_per_gene = 2, verbose = TRUE) {

  counts <- dough$data$counts
  row_data <- dough$data$row_data
  col_data <- dough$data$col_data
  ntc <- dough$data$ntc

  # indicate which guides are NTC
  is_ntc <- row_data$sgRNA %in% ntc$guide

  # filter targeting guides by minimum count across all samples
  keep_counts <- rowSums(counts) >= min_counts

  # filter genes by minimum guide per gene
  gene_table <- table(row_data$gene)
  keep_genes <- row_data$gene %in% names(gene_table[gene_table >= min_guides_per_gene])

  # keep guides that satisfy counts and gene rules, or are NT controls
  keep <- (keep_counts & keep_genes) | is_ntc

  counts_filtered <- counts[keep, , drop = FALSE]
  row_data_filtered <- row_data[keep, , drop = FALSE]

  # reconstruct count matrix so targeting guides come first and then non-targeting guides
  if(!is.null(ntc)) {
    is_ntc_filtered <- row_data_filtered$sgRNA %in% ntc$guide
    counts_filtered <- rbind(counts_filtered[!is_ntc_filtered, , drop = FALSE],
                             counts_filtered[is_ntc_filtered, , drop = FALSE])
    row_data_filtered <- rbind(row_data_filtered[!is_ntc_filtered, , drop = FALSE],
                          row_data_filtered[is_ntc_filtered, , drop = FALSE])
    ntc <- ntc[ntc$guide %in% row_data_filtered$sgRNA, , drop = FALSE]
    ntc$index <- match(ntc$guide, row_data_filtered$sgRNA)
  }

  # save filtered data to dough
  dough$data <- list(
    counts = counts_filtered,
    row_data = row_data_filtered,
    col_data = col_data,
    ntc = ntc
  )

  # print before filtering and after filtering
  if (verbose) {
    message("guides before: ", nrow(row_data),
            "; after: ", nrow(row_data_filtered))
    message("genes before: ", length(unique(row_data$gene)),
            "; after: ", length(unique(row_data_filtered$gene)))
    message("non-targeting guides preserved: ", ifelse(!is.null(nrow(ntc)), nrow(ntc), 'no non-targeting guides were provided'))
  }

  return(dough)
}


