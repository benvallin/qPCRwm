#' Import melt curve data from QuantStudio qPCR software
#'
#'import_melt() tidies the "Melt Curve Raw Data" sheet of an .xls file produced by QuantStudio Design & Analysis Software v1.5.2.
#'
#' @param file character vector of length 1. Path to the .xls file containing the qPCR results to import.
#' @param from character vector of length 1. Indicates which qPCR software the input file originates from. Currently only supports "QuantStudio".
#'
#' @return A tibble with 7 columns: "well_position", "sample_id", "tar_nm", "reading", "temperature", "fluorescence", "derivative".
#'
#' @export
#'
#' @examples
#'
import_melt <- function(file, from = "QuantStudio") {

  # Stop execution in case of invalid input

  if (!requireNamespace("readxl", quietly = TRUE)) {
    warning("The readxl package must be installed to use this function.",
            call. = FALSE)
    return(NULL)
  }

  if (!file.exists(file)) {
    stop("File \"", file, "\" does not exist.",
         call. = FALSE)
  }

  if (length(from) != 1L || !from %in% c("QuantStudio")) {
    stop("from value must be \"QuantStudio\".",
         call. = FALSE)
  }


  invalid_file <- vapply(X = c("Melt Curve Raw Data", "Results"),
                         FUN = function(x) {
                           unique(class(try(expr = match(x = "Well",
                                                         table = suppressMessages(expr = readxl::read_excel(path = file,
                                                                                                            sheet = x))[[1]]),
                                            silent = TRUE)) == "try-error") },
                         FUN.VALUE = vector(mode = "logical", length = 1L))

  if (any(invalid_file)) {
    stop("The input file is not valid.",
         call. = FALSE)
  }

  skip <- vapply(X = c("Melt Curve Raw Data", "Results"),
                 FUN = function(x) { match(x = "Well",
                                           table = suppressMessages(expr = readxl::read_excel(path = file,
                                                                                              sheet = x))[[1]]) },
                 FUN.VALUE = vector(mode = "integer", length = 1L))

  if (any(is.na(skip))) {
    stop("The \"Results\" and/or \"Melt Curve Raw Data\" sheet(s) in the input file do(es) not contain the required \"Well\" column.",
         call. = FALSE)
  }

  # Read qPCR results file

  melt_data <- suppressMessages(expr = readxl::read_excel(path = file,
                                                          sheet = "Melt Curve Raw Data",
                                                          skip = skip[["Melt Curve Raw Data"]]))

  sample_data <- suppressMessages(expr = readxl::read_excel(path = file,
                                                            sheet = "Results",
                                                            skip = skip[["Results"]]))

  # Tidy QuantStudio input

  if (!all(c("Well Position", "Reading", "Temperature", "Fluorescence", "Derivative", "Target Name") %in% colnames(melt_data))) {

    stop("The \"Melt Curve Raw Data\" sheet in the input file does not contain one or more required column(s).",
         call. = FALSE)
  }

  if (!all(c("Well Position", "Sample Name", "Target Name") %in% colnames(sample_data))) {
    stop("The \"Results\" sheet in the input file does not contain one or more required column(s).",
         call. = FALSE)
  }

  melt_data <- stats::setNames(object = melt_data[, c("Well Position", "Reading", "Temperature", "Fluorescence", "Derivative", "Target Name")],
                               nm = c("well_position", "reading", "temperature", "fluorescence", "derivative", "tar_nm"))

  sample_data <- stats::setNames(object = sample_data[, c("Well Position", "Sample Name", "Target Name")],
                                 nm = c("well_position", "sample_id", "tar_nm"))



  # Merge melt curve and sample ID data and tidy final output

  data <- merge(melt_data, sample_data, all = TRUE)

  data <- data[!is.na(data$sample_id), c("well_position", "sample_id", "tar_nm", "reading", "temperature", "fluorescence", "derivative")]

  data <- data[order(data$sample_id, data$tar_nm, data$reading), , drop = FALSE]

  return(data)

}
