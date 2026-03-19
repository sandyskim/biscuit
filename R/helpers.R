#' Calculate size factors via median-of-ratios (ref. DESeq2)
#'
#' @param dough dough object with $data$counts
#' @return size factors
#' @export
compute_size_factors <- function(dough) {
  counts <- dough$data$counts
  log_gm <- rowMeans(log(counts))

  sf <- apply(counts, 2, function(cnts) {
    exp(median((log(cnts) - log_gm)[is.finite(log_gm) & cnts > 0]))
  })

  return(sf)
}

#' Normalize counts using the median-of-ratios method (ref. DESeq2)
#'
#' @param dough dough (or biscuit) object with $data$counts
#' @return matrix of normalized read counts
#' @export
normalize_counts <- function(dough) {
  norm_counts <- sweep(dough$data$counts, 2, compute_size_factors(dough), "/")

  return(norm_counts)
}


#' Compute local false discovery rate (lfdr) and local false sign rate (lfsr) for mu/beta1, adjusted by estimated null distribution
#'
#' @param samples numeric vector of posterior samples
#' @return numeric between 0 and 1
#' @export
compute_lfsr_from_csv <- function(csv_files, par_pattern, stan_variables, mu_ntc_col = NULL, tau_ntc_col = NULL) {
  par_cols <- grep(par_pattern, stan_variables, value = TRUE)
  n_pos_lfdr <- numeric(length(par_cols))
  n_neg_lfdr <- numeric(length(par_cols))
  n_pos_lfsr <- numeric(length(par_cols))
  n_neg_lfsr <- numeric(length(par_cols))
  n_two_sided <- numeric(length(par_cols))
  n_tot <- 0L

  for (f in csv_files) {
    # strip comment lines before passing to fread
    lines <- readLines(f)
    data_lines <- lines[!startsWith(lines, "#")]
    tmp <- tempfile(fileext = ".csv")
    writeLines(data_lines, tmp)
    on.exit(unlink(tmp), add = TRUE)

    select_cols <- if (!is.null(mu_ntc_col) | !is.null(tau_ntc_col)) c(par_cols, mu_ntc_col, tau_ntc_col) else par_cols
    chunk <- data.table::fread(tmp, select = select_cols)
    # drop warmup draws, keeping only sampling iterations
    if ("warmup__" %in% colnames(chunk)) {
      chunk <- chunk[warmup__ == 0]
    }
    mu_draws <- as.matrix(chunk[, par_cols, with = FALSE])
    delta <- if (!is.null(mu_ntc_col)) {
      sweep(mu_draws, 1, chunk[[mu_ntc_col]], "-")
    } else {
      mu_draws
    }

    threshold <- if (!is.null(tau_ntc_col)) chunk[[tau_ntc_col]] else 0
    n_pos_lfdr <- n_pos_lfdr + colSums(delta > threshold)
    n_neg_lfdr <- n_neg_lfdr + colSums(delta < -threshold)
    n_two_sided <- n_two_sided + colSums(abs(delta) < threshold)
    n_pos_lfsr <- n_pos_lfsr + colSums(delta > 0)
    n_neg_lfsr <- n_neg_lfsr + colSums(delta < 0)
    n_tot <- n_tot + nrow(chunk)
  }

  data.frame(
    variable = gsub("^(\\w+)\\.(\\d+)$", "\\1[\\2]", par_cols),
    lfdr = n_two_sided / n_tot,
    lfdr.pos = n_neg_lfdr / n_tot,
    lfdr.neg = n_pos_lfdr / n_tot,
    lfsr = pmin(n_neg_lfsr, n_pos_lfsr) / n_tot,
    lfsr.pos = n_pos_lfsr / n_tot,
    lfsr.neg = n_neg_lfsr / n_tot
  )
}
