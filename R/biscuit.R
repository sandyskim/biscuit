#' generate model input data, intended for internal use
#'
#' @param dough dough object, with data stored in $data
#' @param pseudocount logical indicating whether a pseudocount should be added to the count matrix
#' @return list of input data used to fit the biscuit model
#' @export
generate_biscuit_input <- function(dough, pseudocount=TRUE) {
  # extract data from dough
  counts   <- dough$data$counts
  row_data <- dough$data$row_data
  col_data <- dough$data$col_data
  ntc <- dough$ntc

  # add pseudocount
  if (pseudocount) {
    counts <- counts + 1L
    message("added pseudocount of 1 to counts matrix.")
  }

  # estimate size factors using DESeq2
  sf <- estimateSizeFactorsForMatrix(counts)

  # normalize counts
  norm_counts <- t(t(counts) / sf)

  # code sample conditions to 0 as control, 1 as treatment
  design <- as.integer(factor(col_data$design)) - 1L

  # indicate ntcs
  is_ntc <- as.integer(row_data$sgRNA %in% ntc$guide)

  # code guide to gene mapping as integers for a numeric mapping
  gene_ids <- as.integer(factor(row_data$gene, levels = unique(row_data$gene)))
  n_genes  <- length(unique(gene_ids)) - as.integer(sum(is_ntc) > 0)

  # put together model data for biscuit
  model_data <- list(
    n_samples     = length(design),
    n_guides      = nrow(counts),
    n_genes       = n_genes,
    guide_to_gene = gene_ids,
    sf            = sf,
    mu_g          = log(rowMeans(norm_counts) + 1e-6),
    is_ntc        = is_ntc,
    x             = design,
    y             = counts
  )

  return(model_data)
}

#' runs biscuit
#'
#' @param dough dough object, with data stored in $data
#' @param output_dir directory to save biscuit output files
#' @param save_samples logical to indicate
#' @param n_parallel_chains integer indicating the number chains to run in parallel
#' @param seed integer indicating the seed for reproducibility
#' @param pseudocount logical indicating whether a pseudocount should be added to the count matrix
#' @return biscuit object, with $data and $fit
#' @export
fit_biscuit <- function(dough, output_dir, save_samples=TRUE, n_parallel_chains=4, seed=13, pseudocount=TRUE) {
  # create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir)
  }

  # generate model data input
  model_data <- knead_dough(dough, pseudocount)

  # pick Stan file depending on provided non-targeting controls
  if (is.null(dough$data$ntc)) {
    message("no non-targeting controls detected, using non-targeting control free model.")
    stan_file <- system.file("stan", "crispr_screen_no_ntcs.stan", package = "biscuit")
  } else {
    stan_file <- system.file("stan", "crispr_screen.stan", package = "biscuit")
  }

  # compile stan model
  mod <- cmdstan_model(stan_file)

  # generate reproducible per-chain seeds
  chain_seeds <- seed + seq_len(n_chains) - 1

  # save output into log
  sink(file.path(output_dir, "biscuit.log"))
  on.exit(sink(), add = TRUE)

  # sample
  fit <- mod$sample(
    data = model_data,
    parallel_chains = n_parallel_chains,
    seed = chain_seeds,
    show_exceptions = FALSE,
    show_messages = TRUE
  )

  # create biscuit object
  biscuit <- list(
    data = dough$data,
    fit = list(
      posterior = fit$draws(format='df'),
      diagnostics = fit$diagnostic_summary(),
      runtime = fit$time()
    )
  )
  class(biscuit) <- "biscuit"

  dir.create(file.path(output_dir, "biscuit_output/"), recursive = TRUE, showWarnings = FALSE)

  # save samples if indicated
  if (save_samples) {
    saveRDS(biscuit$fit$posterior, file.path(output_dir, "biscuit_output/posterior.Rdata"))
  }

  return(biscuit)
}

#' generate model input data, intended for internal use (wrapper for generate_biscuit_input)
#'
#' @param dough dough object, with data stored in $data
#' @param pseudocount logical indicating whether a pseudocount should be added to the count matrix
#' @return list of input data used to fit the biscuit model
#' @export
knead_dough <- function(dough, pseudocount=TRUE) {
  model_data <- generate_biscuit_input (dough, pseudocount)
  return(model_data)
}

#' run biscuit (wrapper for run_biscuit)
#'
#' @param dough dough object, with data stored in $data
#' @param output_dir directory to save biscuit output files
#' @param save_samples logical to indicate
#' @param n_parallel_chains integer indicating the number chains to run in parallel
#' @param seed integer indicating the seed for reproducibility
#' @param pseudocount logical indicating whether a pseudocount should be added to the count matrix
#' @return biscuit object, with $data and $fit
#' @export
bake_biscuit <- function(dough, output_dir, save_samples=TRUE, n_parallel_chains=4, seed=13, pseudocount=TRUE) {
  biscuit <- fit_biscuit(dough, output_dir, save_samples, n_parallel_chains, seed, pseudocount)
  return(biscuit)
}

