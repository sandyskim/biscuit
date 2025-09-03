#' compute local false sign rate (lfsr)
#'
#' @param samples numeric vector of posterior samples
#' @return numeric between 0 and 0.5
#' @export
compute_lfsr <- function(samples) {
  p_pos <- mean(samples >= 0)
  p_neg <- mean(samples <= 0)
  lfsr <- min(p_pos, p_neg)

  return(lfsr)
}

#' extract posterior draws for one parameter
#'
#' @param biscuit fitted biscuit object
#' @param pars single parameter base name (e.g. "mu")
#' @return tibble with columns: draw, parameter index, value
#' @export
extract_parameters <- function(biscuit, pars) {
  draws <- as.data.frame(biscuit$fit$posterior)

  # select only columns for this parameter
  matching_cols <- grep(
    paste0("^", pars, "\\["),
    colnames(draws),
    value = TRUE
  )
  if (length(matching_cols) == 0) return(tibble())

  parameters <- as_tibble(draws[, matching_cols, drop = FALSE]) %>%
    mutate(draw = row_number()) %>%
    pivot_longer(
      cols = all_of(matching_cols),
      names_to = "parameter",
      values_to = "value"
    ) %>%
    mutate(
      index = as.integer(str_extract(parameter, "(?<=\\[)\\d+(?=\\])"))
    ) %>%
    select(draw, index, value)

  return(parameters)
}


#' summarize posterior draws
#'
#' @param biscuit fitted biscuit object
#' @param pars list parameters to summarize (e.g. mu, beta1)
#' @return biscuit object with $results
#' @export
summarize_params <- function(biscuit,
                             pars = c("beta1", "mu", "beta0", "phi", "gamma")) {
  biscuit$results <- list()

  for (par in pars) {
    draws <- extract_parameters(biscuit, pars = par)
    if (nrow(draws) == 0) next

    summ <- draws %>%
      group_by(index) %>%
      summarize(
        mean   = mean(value),
        median = median(value),
        sd     = sd(value),
        q2.5   = quantile(value, 0.025),
        q97.5  = quantile(value, 0.975),
        lfsr   = if (par %in% c("beta1", "mu")) compute_lfsr(value) else NA_real_,
        .groups = "drop"
      ) %>%
      arrange(index)

    if (par %in% c("beta1", "beta0", "phi")) {
      summ <- summ %>%
        mutate(
          sgRNA = biscuit$data$row_data$sgRNA[index],
          gene  = biscuit$data$row_data$gene[index]
        )
    } else if (par == "mu") {
      genes <- unique(biscuit$data$row_data$gene)
      summ <- summ %>%
        mutate(gene = genes[index])
    } else if (par == "gamma") {
      summ <- summ %>%
        mutate(sample = biscuit$data$col_data$sample[index])
    }

    biscuit$results[[par]] <- summ
  }

  return(biscuit)
}
