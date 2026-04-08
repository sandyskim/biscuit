#' Generate model input data for the biscuit Stan model
#'
#' @param dough dough object with \code{$data}
#' @param pseudocount logical; whether to add a pseudocount of 1 to the count matrix
#' @return named list of input data for the biscuit Stan model
#' @export
generate_biscuit_input <- function(dough, pseudocount = TRUE) {
  counts <- dough$data$counts
  row_data <- dough$data$row_data
  col_data <- dough$data$col_data
  controls <- dough$data$controls

  if (pseudocount) {
    counts <- counts + 1L
    message("Added pseudocount of 1 to counts matrix.")
  }

  sf <- compute_size_factors(dough)
  norm_counts <- sweep(counts, 2, sf, "/")
  design <- as.integer(col_data$design) - 1L
  is_ntc <- as.integer(row_data$sgRNA %in% controls$guide)

  unique_genes <- unique(row_data$gene[!as.logical(is_ntc)])
  gene_ids <- as.integer(factor(row_data$gene, levels = unique(row_data$gene)))
  n_genes <- length(unique_genes)

  control_cols <- which(col_data$design == "control")
  mean_control <- rowMeans(norm_counts[, control_cols, drop = FALSE])
  mu_g <- as.numeric(log(rowMeans(norm_counts + 1)))
  beta0_hat <- if (length(control_cols) == 1) mu_g else as.numeric(log(mean_control + 1))

  list(
    n_samples = length(design),
    n_guides = nrow(counts),
    n_genes = n_genes,
    guide_to_gene = gene_ids,
    y = counts,
    sf = sf,
    beta0_hat = beta0_hat,
    mu_g = mu_g,
    is_ntc = is_ntc,
    x = design
  )
}


#' Fit the biscuit Stan model
#'
#' @param dough dough object with data stored in \code{$data}
#' @param output_dir directory to save biscuit output files
#' @param filter logical; whether to filter counts before fitting
#' @param save_samples logical; whether to save posterior draws
#' @param n_parallel_chains integer; number of chains to run in parallel
#' @param seed integer; seed for reproducibility
#' @param pseudocount logical; whether to add a pseudocount to the count matrix
#' @return biscuit object with \code{$data} and \code{$fit}
#' @export
fit_biscuit <- function(dough,
                        output_dir,
                        filter = TRUE,
                        save_samples = TRUE,
                        n_parallel_chains = 4,
                        seed = 13,
                        pseudocount = TRUE) {
  if (!dir.exists(output_dir)) dir.create(output_dir)

  if (filter) {
    message("Filtering counts...")
    dough <- trim_dough(dough)
  }

  model_data <- generate_biscuit_input(dough, pseudocount)

  stan_file <- if (is.null(dough$data$controls)) {
    message("No non-targeting controls detected, using NTC-free model.")
    system.file("stan", "crispr_screen_no_ntcs.stan", package = "biscuit")
  } else {
    system.file("stan", "crispr_screen.stan", package = "biscuit")
  }

  mod <- cmdstan_model(stan_file)

  sink(file.path(output_dir, "biscuit.log"), split = TRUE)
  on.exit(sink(), add = TRUE)

  fit <- mod$sample(
    data = model_data,
    parallel_chains = n_parallel_chains,
    seed = seed,
    show_exceptions = FALSE,
    show_messages = TRUE
  )

  # summarize and add lfdr and lfsr
  fit_summary <- fit$summary()
  rownames(fit_summary) <- fit_summary$variable

  csv_files <- fit$output_files()
  stan_variables <- readLines(csv_files[1]) |>
    (\(lines) lines[!startsWith(lines, "#")])() |>
    (\(lines) lines[1])() |>
    strsplit(",") |>
    unlist()
  lfsr_all <- rbind(
    compute_lfsr_from_csv(csv_files, "^mu\\.", stan_variables, mu_ntc_col = "mu_ntc", tau_ntc_col = "tau_ntc"),
    compute_lfsr_from_csv(csv_files, "^beta1\\.", stan_variables)
  )
  fit_summary <- merge(fit_summary, lfsr_all, by = "variable", all.x = TRUE)

  results = summarize_parameters(dough, fit_summary)

  # save outputs
  output_path <- file.path(output_dir, "biscuit_output")
  dir.create(output_path, recursive = TRUE, showWarnings = FALSE)
  write.csv(fit_summary, file.path(output_path, "fit_summary.csv"))

  posterior_path <- NULL
  if (save_samples) {
    posterior_path <- file.path(output_path, "_posterior.csv")
    con <- file(posterior_path, open = "wt")
    on.exit(close(con), add = TRUE)

    writeLines(grep("^#", readLines(csv_files[1]), invert = TRUE, value = TRUE)[1], con)
    for (f in csv_files) {
      lines <- readLines(f)
      data_lines <- grep("^#", lines, invert = TRUE, value = TRUE)[-1]
      writeLines(data_lines, con)
    }
  }

  list(
      data = dough$data,
      fit  = list(
        diagnostics = fit$diagnostic_summary(),
        runtime = fit$time(),
        model = fit$code(),
        model_data = model_data,
        posterior_path = posterior_path
      ),
      results = results
  )
}
