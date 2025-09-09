#' extract posterior draws for one parameter
#'
#' @param biscuit fitted biscuit object
#' @param pars parameter name (e.g. "mu")
#' @return data frame with columns: draw, parameter index, value
#' @export
extract_parameters <- function(biscuit, pars) {
  draws <- as.data.frame(biscuit$fit$posterior)

  # select only columns for this parameter
  matching_cols <- grep(paste0("^", pars, "\\["),
                        colnames(draws),
                        value = TRUE)
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
summarize_parameters <- function(biscuit, pars = c("mu", "beta1", "beta0", "phi", "gamma")) {
    biscuit$results <- lapply(pars, function(par) {
      draws <- extract_parameters(biscuit, par)
      if (nrow(draws) == 0)
        return(NULL)

      summ <- draws %>%
        group_by(index) %>%
        summarize(
          mean   = mean(value),
          median = median(value),
          sd     = sd(value),
          q2.5   = quantile(value, 0.025),
          q97.5  = quantile(value, 0.975),
          lfsr   = if (par %in% c("beta1", "mu"))
            compute_lfsr(value)
          else
            NA_real_,
          .groups = "drop"
        ) %>%
        arrange(index)

      # add identifiers
      if (par %in% c("beta1", "beta0", "phi")) {
        summ <- summ %>%
          mutate(
            sgRNA = biscuit$data$row_data$sgRNA[index],
            gene  = biscuit$data$row_data$gene[index]
          )
      } else if (par == "mu") {
        genes <- unique(biscuit$data$row_data$gene)
        summ <- summ %>% mutate(gene = genes[index])
      } else if (par == "gamma") {
        summ <-
          summ %>% mutate(sample = biscuit$data$col_data$sample[index])
      }

      summ
    })

    names(biscuit$results) <- pars

    return(biscuit)
  }
