#' calculate size factors via median-of-ratios (ref. DESeq2)
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


#' normalize counts using the median-of-ratios method (ref. DESeq2)
#'
#' @param dough dough (or biscuit) object with $data$counts
#' @return matrix of normalized read counts
#' @export
normalize_counts <- function(dough) {
  counts <- dough$data$counts
  sf <- compute_size_factors(dough)
  norm_counts <- t(t(counts) / sf)

  return(norm_counts)
}

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
