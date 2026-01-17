#' extract posterior draws for one parameter
#'
#' @param biscuit fitted biscuit object
#' @param pars parameter name (e.g. "mu")
#' @return data frame with columns: draw, parameter index, value
#' @export
extract_parameters <- function(biscuit, pars) {
  draws <- tibble::as_tibble(biscuit$fit$posterior)

  matching_cols <- grep(paste0("^", pars, "(\\[|$)"), colnames(draws), value = TRUE)
  if (length(matching_cols) == 0)
    return(data.frame())

  parameters <- draws %>%
    select(all_of(matching_cols)) %>%
    mutate(draw = row_number()) %>%
    pivot_longer(cols = -draw,
                 names_to = "parameter",
                 values_to = "value") %>%
    mutate(index = as.integer(str_extract(parameter, "(?<=\\[)\\d+(?=\\])"))) %>%
    select(draw, index, value)

  return(parameters)
}

#' summarize posterior draws
#'
#' @param biscuit biscuit object with $fit
#' @param pars list parameter names to summarize (e.g. c('mu', 'beta1'))
#' @return biscuit object with $results
#' @export
summarize_parameters <- function(biscuit, pars = c("mu", "beta1", "beta0", "phi", "eff", "mu_ntc", "tau_ntc")) {

  posterior <- as.data.frame(biscuit$fit$posterior)
  row_data <- biscuit$data$row_data
  col_data <- biscuit$data$col_data
  genes_unique <- unique(row_data$gene)

  mu_ntc_cols <- grep("^mu_ntc", colnames(posterior))
  tau_cols    <- grep("^tau_ntc", colnames(posterior))

  mu_ntc_draws <- if (length(mu_ntc_cols)) posterior[[mu_ntc_cols[1]]] else rep(0, nrow(posterior))
  tau_draws    <- if (length(tau_cols))    posterior[[tau_cols[1]]]    else rep(0, nrow(posterior))

  biscuit$results <- lapply(pars, function(par) {
    draws <- extract_parameters(biscuit, par)
    if (nrow(draws) == 0) return(NULL)
    summ <- draws %>%
      group_by(index) %>%
      summarize(
        mean   = mean(value),
        median = median(value),
        sd     = sd(value),
        q2.5   = quantile(value, 0.025),
        q97.5  = quantile(value, 0.975),
        .groups = "drop"
      ) %>%
      arrange(index)
    if (par %in% c("beta1", "beta0", "phi")) {
      summ$sgRNA <- row_data$sgRNA[summ$index]
      summ$gene  <- row_data$gene[summ$index]
    }
    if (par == "beta1") {
      beta1_cols     <- grep("^beta1\\[", colnames(posterior))
      if (!is.null(mu_ntc_draws) && !is.null(tau_draws)) {
        beta1_draws_matrix <- posterior[, beta1_cols, drop = FALSE]
        lfdr <- vapply(seq_len(ncol(beta1_draws_matrix)), function(j) {
          compute_rope_lfdr(
            mu_g   = beta1_draws_matrix[[j]],
            mu_ntc = mu_ntc_draws,
            tau    = tau_draws
          )
        }, FUN.VALUE = numeric(1))
        summ$lfdr <- lfdr
      } else {
        summ$lfdr <- NA_real_
      }
    }
    if (par == "mu") {
      summ$gene <- genes_unique[summ$index]
      mu_cols     <- grep("^mu\\[", colnames(posterior))
      mu_draws_matrix <- posterior[, mu_cols, drop = FALSE]
      lfsr <- vapply(seq_len(ncol(mu_draws_matrix)), function(j) {
        compute_lfsr(
          mu_g   = mu_draws_matrix[[j]],
          mu_ntc = mu_ntc_draws,
        )
      }, FUN.VALUE = numeric(1))
      summ$lfsr <- lfsr
      lfsr.neg <- vapply(seq_len(ncol(mu_draws_matrix)), function(j) {
        compute_lfsr(
          mu_g   = mu_draws_matrix[[j]],
          mu_ntc = mu_ntc_draws,
          mode = "neg"
        )
      }, FUN.VALUE = numeric(1))
      summ$lfsr.neg <- lfsr.neg
      lfsr.pos <- vapply(seq_len(ncol(mu_draws_matrix)), function(j) {
        compute_lfsr(
          mu_g   = mu_draws_matrix[[j]],
          mu_ntc = mu_ntc_draws,
          mode = "pos"
        )
      }, FUN.VALUE = numeric(1))
      summ$lfsr.pos <- lfsr.pos
        lfdr <- vapply(seq_len(ncol(mu_draws_matrix)), function(j) {
          compute_rope_lfdr(
            mu_g   = mu_draws_matrix[[j]],
            mu_ntc = mu_ntc_draws,
            tau    = tau_draws
          )
        }, FUN.VALUE = numeric(1))
        summ$lfdr <- lfdr
        lfdr.neg <- vapply(seq_len(ncol(mu_draws_matrix)), function(j) {
          compute_rope_lfdr(
            mu_g   = mu_draws_matrix[[j]],
            mu_ntc = mu_ntc_draws,
            tau    = tau_draws,
            mode = 'neg'
          )
        }, FUN.VALUE = numeric(1))
        summ$lfdr.neg <- lfdr.neg
        lfdr.pos <- vapply(seq_len(ncol(mu_draws_matrix)), function(j) {
          compute_rope_lfdr(
            mu_g   = mu_draws_matrix[[j]],
            mu_ntc = mu_ntc_draws,
            tau    = tau_draws,
            mode = 'pos'
          )
        }, FUN.VALUE = numeric(1))
        summ$lfdr.pos <- lfdr.pos
    }
    summ
  })

  names(biscuit$results) <- pars
  return(biscuit)
}
