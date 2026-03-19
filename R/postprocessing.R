#' Summarize posterior parameters
#'
#' @param dough dough object with \code{$fit}
#' @param pars character vector of parameter names to summarize
#' @return dough object with \code{$results} appended
#' @export
summarize_parameters <- function(dough, fit_summary,
                                 pars = c("mu", "beta1", "beta0", "phi", "eff", "mu_ntc", "tau_ntc")) {
  row_data <- dough$data$row_data
  genes_unique <- unique(row_data$gene)
  fit_summary <- fit_summary

  stats::setNames(lapply(pars, function(par) {
    summ <- fit_summary[grep(paste0("^", par, "\\["), fit_summary$variable), ]
    if (nrow(summ) == 0) {
      summ <- fit_summary[fit_summary$variable == par, ]
    }
    if (nrow(summ) == 0) return(NULL)
    if(nrow(summ) > 1) {
      summ$index <- as.integer(gsub(".*\\[(\\d+)\\]$", "\\1", summ$variable))
    }
    else {
      summ$index <- 1
    }
    summ <- summ[order(summ$index), ]

    stat_cols <- grep("^lfsr|^lfdr", colnames(summ), value = TRUE)
    if (!par %in% c("mu", "beta1") && length(stat_cols)) {
      summ <- summ[, !colnames(summ) %in% stat_cols, drop = FALSE]
    }

    if (par %in% c("beta1", "beta0", "phi")) {
      summ$sgRNA <- row_data$sgRNA[summ$index]
      summ$gene <- row_data$gene[summ$index]
    }

    if (par == "mu") {
      summ$gene <- genes_unique[summ$index]
    }

    summ
  }), pars)
}
