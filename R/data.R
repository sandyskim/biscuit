#' REGNASE-1 CRISPR screen data
#'
#' contains results from a CRISPR-Cas9 screen targeting Regnase-1 in OT-I T cells.
#'
#' @format list with lists
#' \describe{
#'   \item{counts}{n_guides rows, n_samples columns matrix: raw sequencing counts from CRISPR screen}
#'   \item{guide_to_gene}{n_guides rows, 2 columns matrix: names of guides and the gene they target}
#'   \item{sample_design}{n_samples length vector: sample conditions}
#'   \item{ntcs}{n_guides length vector: names of non-targeting control guides}
#' }
#'
#' @source Wei, J., Long, L., Zheng, W., et al. (2019). Targeting Regnase-1 programs long-lived effector T cells for cancer therapy. Nature, 576(7787), 471-476. doi:10.1038/s41586-019-1821-z
"regnase"

