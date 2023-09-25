#' Compute the expression fold changes for every gene of interest
#'
#'compute_fold_change() calculates gene expression fold changes, provided a data.frame, a grouping variable, a calibrator sample and one or more reference gene(s).
#'
#' @param data data.frame / tibble. Must contain the columns group_var, "tar_nm" and ct_var.
#' @param group_var character vector of length 1. Name of the column to use as grouping variable. Default is "sample_id".
#' @param ct_var character vector of length 1. Name of the column containing the Ct values. Default is "mean_id".
#' @param ref_nm character vector. Name(s) of the gene(s) to use as reference gene(s). All ref_nm elements must be in the "tar_nm" column of the input data.frame.
#' @param cal_nm character vector of length 1. Name of the sample to use as calibrator. Must be in the group_var column of the input data.frame.
#' @param method character vector of length 1. Indicates which method to use for computation of gene expression fold changes. Valid values are "ddct" or "Pfaffl".
#' @param use_geomean logical vector of length 1. Indicates whether to use geometric mean instead of arithmetic mean for computing the mean of multiple reference genes Ct values. Ignored if method is "Pfaffl". Default is FALSE.
#' @param efficiencies data.frame / tibble. Table with gene-specific primer efficiencies. Must contain the columns "tar_nm" and "efficiency", and one entry for each gene in the "tar_nm" column of the input data.frame. Ignored if method is "ddct". If missing and method is "Pfaffl", efficiency is set to 2 for every gene.
#'
#' @return A data.frame of variable length: 10 (ddct method), 12 (Pfaffl method and multiple ref_nm), or 14 (Pfaffl method and single ref_nm). Gene expression fold changes are stored in the "fold_change" column.
#'
#' @export
#'
#' @importFrom stats setNames
#'
#' @examples
#'
compute_fold_change <- function(data,
                                group_var = "sample_id",
                                ct_var = "mean_ct",
                                ref_nm,
                                cal_nm,
                                method = "ddct",
                                use_geomean = FALSE,
                                efficiencies = NULL) {

  # Capture user's arguments

  data <- eval(expr = substitute(data), envir = parent.frame())

  group_var <- eval(expr = substitute(group_var), envir = parent.frame())

  ct_var <- eval(expr = substitute(ct_var), envir = parent.frame())

  ref_nm <- eval(expr = substitute(ref_nm), envir = parent.frame())

  cal_nm <- eval(expr = substitute(cal_nm), envir = parent.frame())

  method <- eval(expr = substitute(method), envir = parent.frame())

  use_geomean <- eval(expr = substitute(use_geomean), envir = parent.frame())

  efficiencies <- eval(expr = substitute(efficiencies), envir = parent.frame())

  # Determine if quantification includes a single or multiple reference gene(s)

  multi_ref_nm <- ifelse(length(ref_nm) > 1L, T, F)

  # Stop execution in case of invalid input

  if (length(group_var) != 1 || !is.character(group_var)) {
    stop("\ngroup_var value must be a character vector of length 1.\n")
  }

  if (length(ct_var) != 1 || !is.character(ct_var)) {
    stop("\nct_var value must be a character vector of length 1.\n")
  }

  if (!is.data.frame(data) || !all(c(group_var, "tar_nm", ct_var) %in% colnames(data))) {
    stop("\ndata must be a dataframe containing the columns \"", group_var, "\", \"tar_nm\" and \"", ct_var, "\".\n")
  }

  if (any(is.na(data[[group_var]]))) {
    stop("\nThe group_var column in the input data must contain only non-NA values.\n")
  }

  if (!is.character(data[[group_var]])) {
    stop("\nThe group_var column in the input data must be of class character.\n")
  }

  if (!is.numeric(data[[ct_var]])) {
    stop("\nThe ct_var column in the input data must be of class numeric.\n")
  }

  if (any(is.na(data$tar_nm))) {
    stop("\nThe \"tar_nm\" column in the input data must contain only non-NA values.\n")
  }

  if (!is.character(data$tar_nm)) {
    stop("\nThe \"tar_nm\" column in the input data must be of class character.\n")
  }

  if (!is.character(ref_nm)) {
    stop("\nref_nm value must be a character vector.\n")
  }

  if (!all(ref_nm %in% data$tar_nm)) {
    stop("\nref_nm value(s) must be in the \"tar_nm\" column of the input data.\n")
  }

  if (length(cal_nm) != 1 || !is.character(cal_nm)) {
    stop("\ncal_nm value must be a character vector of length 1.\n")
  }

  if (!cal_nm %in% data[[group_var]]) {
    stop("\ncal_nm value must be in the \"", group_var, "\" column of the input data.\n")
  }

  if (length(method) != 1L || !method %in% c("ddct", "Pfaffl")) {
    stop("\nmethod value must be \"ddct\" or \"Pfaffl\".\n")
  }

  if (!is.null(efficiencies)) {
    if (!is.data.frame(efficiencies) || !all(c("tar_nm", "efficiency") %in% colnames(efficiencies))) {
      stop("\nefficiencies must be a dataframe containing the columns \"tar_nm\" and \"efficiency\".\n")
    }
    if (!all(data$tar_nm %in% efficiencies$tar_nm)) {
      stop("\nThe input efficiencies table does not contain an efficiency value for every tar_nm.\n")
    }
    if (!is.numeric(efficiencies$efficiency)) {
      stop("\nThe \"efficiency\" column of efficiencies must be of class numeric.\n")
    }
  }

  if (!use_geomean %in% c(T, F)) {
    stop("\nuse_geomean value must be a logical TRUE or FALSE.\n")
  }

  # Message in case of ignored argument

  if (method == "ddct" && !is.null(efficiencies)) {
    message("\nThe \"ddct\" method is incompatible with primer-specific efficiencies.\nProvided efficiencies will be ignored.\n")
  }

  if (use_geomean && !multi_ref_nm) {
    message("\nThe use_geomean argument serves only in case of multiple reference genes.\nuse_geomean = TRUE will be ignored.\n")
  }

  if (use_geomean && multi_ref_nm && method == "Pfaffl") {
    message("\nThe use_geomean argument is not used by the \"Pfaffl\" method (does not compute the mean of multiple reference genes Ct values).\nuse_geomean = TRUE will be ignored.\n")
  }

  # Isolate Ct values for target and references genes, respectively
  # In case of multiple reference genes, compute the arithmetic / geometric mean of reference genes Ct values per sample ID

  tar_ct_data <- data[!data$tar_nm %in% ref_nm, c(group_var, "tar_nm", ct_var)]

  indiv_ref_ct_data <- data[data$tar_nm %in% ref_nm, c(group_var, "tar_nm", ct_var)]

  if (multi_ref_nm) {

    ref_ct_data <- split(x = indiv_ref_ct_data, f = indiv_ref_ct_data[[group_var]])

    ref_ct_data <- lapply(X = ref_ct_data,
                          FUN = function(x) {
                            setNames(object = list(x[[group_var]][[1]],
                                                   paste0(ref_nm, collapse = " / "),
                                                   ifelse(use_geomean,
                                                          exp(mean(log(x[[ct_var]]))),
                                                          mean(x[[ct_var]]))),
                                     nm = c(group_var, "tar_nm", ct_var))
                          })

    ref_ct_data <- Reduce(f = function(...) { merge(..., all = T) },
                          x = ref_ct_data)

  } else {

    ref_ct_data <- indiv_ref_ct_data

  }

  ct_data <- merge(tar_ct_data, ref_ct_data, all = T)

  if (method == "ddct") {

    tar_ct_data <- split(x = tar_ct_data, f = tar_ct_data[[group_var]])

    fc_data <- lapply(X = names(tar_ct_data),
                      FUN = function(x) {
                        setNames(object = list(x,
                                               tar_ct_data[[x]]$tar_nm,
                                               paste0(ref_nm, collapse = " / "),
                                               tar_ct_data[[x]][[ct_var]] - ref_ct_data[ref_ct_data[[group_var]] == x, ct_var]),
                                 nm = c(group_var, "tar_nm", "ref_nm", "delta_ct"))
                      })

    tar_ct_data <- Reduce(f = function(...) { merge(..., all = T) },
                          x = tar_ct_data)

    fc_data <- Reduce(f = function(...) { merge(..., all = T) },
                      x = fc_data)

    fc_data$delta_delta_ct <- fc_data$delta_ct - fc_data[fc_data[[group_var]] == cal_nm, "delta_ct"]

    fc_data$two_power_minus_delta_delta_ct <- 2^(-fc_data$delta_delta_ct)

    fc_data$fold_change <- round(x = fc_data$two_power_minus_delta_delta_ct, digits = 12)

    fc_data$log2_fold_change <- log(x = fc_data$fold_change, base = 2)

    fc_data <- merge(fc_data,
                     setNames(object = ct_data[ct_data$tar_nm %in% unique(fc_data$tar_nm),],
                              nm = c(group_var, "tar_nm", paste0("tar_", ct_var))),
                     all = T)

    fc_data <- merge(fc_data,
                     setNames(object = ct_data[ct_data$tar_nm %in% unique(fc_data$ref_nm),],
                              nm = c(group_var, "ref_nm", paste0("ref_", ct_var))),
                     all = T)

    fc_data <- fc_data[, c(group_var, "tar_nm", "ref_nm", paste0("tar_", ct_var), paste0("ref_", ct_var),
                           "delta_ct", "delta_delta_ct", "two_power_minus_delta_delta_ct", "fold_change", "log2_fold_change")]

  }

  if (method == "Pfaffl") {

    if (is.null(efficiencies)) {

      efficiencies <- data.frame(tar_nm = unique(data$tar_nm),
                                 efficiency = rep(2, times = length(unique(data$tar_nm))))

    }

    tar_rq_data <- tar_ct_data

    tar_rq_data$tar_delta_ct <- tar_rq_data[tar_rq_data[[group_var]] == cal_nm, ct_var] - tar_rq_data[[ct_var]]

    tar_rq_data <- merge(tar_rq_data, efficiencies, all.x = T)

    tar_rq_data$eff_power_tar_delta_ct <- tar_rq_data$efficiency ^ tar_rq_data$tar_delta_ct

    tar_rq_data <- setNames(object = tar_rq_data[, c(group_var, "tar_nm", ct_var, "tar_delta_ct", "efficiency", "eff_power_tar_delta_ct")],
                            nm = c(group_var, "tar_nm", paste0("tar_", ct_var), "tar_delta_ct", "tar_efficiency", "eff_power_tar_delta_ct"))

    if (multi_ref_nm) {

      indiv_ref_rq_data <- indiv_ref_ct_data

      indiv_ref_rq_data$ref_delta_ct <- indiv_ref_rq_data[indiv_ref_rq_data[[group_var]] == cal_nm, ct_var] - indiv_ref_rq_data[[ct_var]]

      indiv_ref_rq_data <- merge(indiv_ref_rq_data, efficiencies, all.x = T)

      indiv_ref_rq_data$eff_power_ref_delta_ct <- indiv_ref_rq_data$efficiency ^ indiv_ref_rq_data$ref_delta_ct

      ref_rq_data <- split(x = indiv_ref_rq_data, f = indiv_ref_rq_data[[group_var]])

      ref_rq_data <- lapply(X = names(ref_rq_data),
                            FUN = function(x) {
                              setNames(object = list(x,
                                                     paste0(ref_nm, collapse = " / "),
                                                     exp(mean(log(ref_rq_data[[x]]$eff_power_ref_delta_ct)))),
                                       nm = c(group_var, "ref_nm", "geomean_eff_power_ref_delta_ct"))
                            })

      ref_rq_data <- Reduce(f = function(...) { merge(..., all = T) },
                            x = ref_rq_data)

      indiv_ref_rq_data <- split(x = indiv_ref_rq_data,
                                 f = indiv_ref_rq_data[[group_var]])

      indiv_ref_rq_data <- lapply(X = indiv_ref_rq_data,
                                  FUN = function(x) {
                                    setNames(object = x[, c("tar_nm", ct_var, "ref_delta_ct", "efficiency", "eff_power_ref_delta_ct")],
                                             nm = c("ref_nm", paste0("ref_", ct_var), "ref_delta_ct", "ref_efficiency", "eff_power_ref_delta_ct"))
                                  })

      indiv_ref_rq_data <- cbind(indiv_ref_rq_data)

      indiv_ref_rq_data <- setNames(object = data.frame(group_var = rownames(indiv_ref_rq_data),
                                                        ref_data = indiv_ref_rq_data),
                                    nm = c(group_var, "ref_data"))

      ref_rq_data <- merge(ref_rq_data, indiv_ref_rq_data, all = T)

      fc_data <- merge(tar_rq_data, ref_rq_data, all = T)

      fc_data$eff_power_tar_delta_ct_over_geomean_eff_power_ref_delta_ct <- fc_data$eff_power_tar_delta_ct / fc_data$geomean_eff_power_ref_delta_ct

      fc_data$fold_change <- round(x = fc_data$eff_power_tar_delta_ct_over_geomean_eff_power_ref_delta_ct, digits = 12)

      fc_data$log2_fold_change <- log(x = fc_data$fold_change, base = 2)

      fc_data <- fc_data[, c(group_var, "tar_nm", "ref_nm",
                             paste0("tar_", ct_var), "tar_delta_ct", "tar_efficiency", "eff_power_tar_delta_ct",
                             "ref_data", "geomean_eff_power_ref_delta_ct",
                             "eff_power_tar_delta_ct_over_geomean_eff_power_ref_delta_ct", "fold_change", "log2_fold_change")]

    } else {

      ref_rq_data <- ref_ct_data

      ref_rq_data$ref_delta_ct <- ref_rq_data[ref_rq_data[[group_var]] == cal_nm, ct_var] - ref_rq_data[[ct_var]]

      ref_rq_data <- merge(ref_rq_data, efficiencies, all.x = T)

      ref_rq_data$eff_power_ref_delta_ct <- ref_rq_data$efficiency ^ ref_rq_data$ref_delta_ct

      ref_rq_data <- setNames(object = ref_rq_data[, c(group_var, "tar_nm", ct_var, "ref_delta_ct", "efficiency", "eff_power_ref_delta_ct")],
                              nm = c(group_var, "ref_nm", paste0("ref_", ct_var), "ref_delta_ct", "ref_efficiency", "eff_power_ref_delta_ct"))

      fc_data <- merge(tar_rq_data, ref_rq_data, all = T)

      fc_data$eff_power_tar_delta_ct_over_eff_power_ref_delta_ct <- fc_data$eff_power_tar_delta_ct / fc_data$eff_power_ref_delta_ct

      fc_data$fold_change <- round(x = fc_data$eff_power_tar_delta_ct_over_eff_power_ref_delta_ct, digits = 12)

      fc_data$log2_fold_change <- log(x = fc_data$fold_change, base = 2)

      fc_data <- fc_data[, c(group_var, "tar_nm", "ref_nm",
                             paste0("tar_", ct_var), "tar_delta_ct", "tar_efficiency", "eff_power_tar_delta_ct",
                             paste0("ref_", ct_var), "ref_delta_ct", "ref_efficiency", "eff_power_ref_delta_ct",
                             "eff_power_tar_delta_ct_over_eff_power_ref_delta_ct", "fold_change", "log2_fold_change")]

    }

    fc_data <- fc_data[order(fc_data[[group_var]], fc_data$tar_nm), , drop = F]

  }

  return(fc_data)

}
