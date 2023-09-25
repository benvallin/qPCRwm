#' Compute the mean Ct value for each level of a grouping variable
#'
#'compute_mean_ct() calculates the arithmetic or geometric mean of Ct values, provided a data.frame and a grouping variable.
#'
#' @param data data.frame / tibble. Must contain the columns group_var, "tar_nm" and "ct".
#' @param group_var character vector of length 1. Name of the column to use as grouping variable. Default is "sample_id".
#' @param use_geomean logical vector of length 1. Indicates whether to compute geometric mean instead of arithmetic mean. Default is FALSE.
#'
#' @return A data.frame with 7 columns: group_var, "tar_nm", "raw_data", "n_rep", "n_filt_rep", "mean_ct" and "sd_ct".
#'
#' @export
#'
#' @importFrom stats setNames
#' @importFrom stats na.omit
#' @importFrom stats sd
#'
#' @examples
#'
compute_mean_ct <- function(data, group_var = "sample_id", use_geomean = FALSE) {

  # Capture user's arguments

  data <- eval(expr = substitute(data), envir = parent.frame())

  group_var <- eval(expr = substitute(group_var), envir = parent.frame())

  use_geomean <- eval(expr = substitute(use_geomean), envir = parent.frame())

  # Stop execution in case of invalid input

  if (length(group_var) != 1 || !is.character(group_var)) {
    stop("\ngroup_var value must be a character vector of length 1.\n")
  }

  if (!is.data.frame(data) || !all(c(group_var, "tar_nm", "ct") %in% colnames(data))) {
    stop("\ndata must be a dataframe containing the columns \"", group_var, "\", \"tar_nm\" and \"ct\".\n")
  }

  if (any(is.na(data[[group_var]]))) {
    stop("\nThe group_var column in the input data must contain only non-NA values.\n")
  }

  if (!use_geomean %in% c(T, F)) {
    stop("\nuse_geomean value must be a logical TRUE or FALSE.\n")
  }

  # Construct group_var and tar_nm values for final output

  data_sample_ids <- rep(levels(as.factor(data[[group_var]])),
                         times = length(unique(data$tar_nm)))

  data_tar_nms <- rep(levels(as.factor(data$tar_nm)),
                      each = length(unique(data[[group_var]])))

  # Construct ct_data
  # => Compute arithmetic / geometric mean and SD of Ct values across technical replicates

  ct_data <- data[, c(group_var, "tar_nm", "ct")]

  ct_data <- split(x = ct_data, f = list(ct_data[[group_var]], ct_data$tar_nm))

  ct_data <- mapply(FUN = function(x, y, z) {
    setNames(object = list(x,
                           y,
                           length(z$ct),
                           length(na.omit(z$ct)),
                           ifelse(use_geomean,
                                  ifelse(is.nan(exp(mean(log(z$ct), na.rm = T))),
                                         NA_real_,
                                         exp(mean(log(z$ct), na.rm = T))),
                                  ifelse(is.nan(mean(z$ct, na.rm = T)),
                                         NA_real_,
                                         mean(z$ct, na.rm = T))),
                           sd(z$ct, na.rm = T)),
             nm = c(group_var, "tar_nm", "n_rep", "n_filt_rep", "mean_ct", "sd_ct"))
  },
  data_sample_ids, data_tar_nms, ct_data,
  SIMPLIFY = F)

  ct_data <- Reduce(f = function(...) { merge(..., all = T) },
                    x = ct_data)

  # Construct raw_data

  raw_data <- cbind(split(x = data,
                          f = list(data[[group_var]], data$tar_nm)))

  raw_data <- setNames(object = data.frame(group_var = data_sample_ids,
                                           tar_nm = data_tar_nms,
                                           raw_data = raw_data),
                       nm = c(group_var, "tar_nm", "raw_data"))

  # Merge ct_data and raw_data into data

  data <- merge(ct_data, raw_data, all = T)

  data <- data[, c(group_var, "tar_nm", "raw_data", "n_rep", "n_filt_rep", "mean_ct", "sd_ct")]

  data <- data[order(data[[group_var]], data$tar_nm), , drop = F]

  return(data)

}
