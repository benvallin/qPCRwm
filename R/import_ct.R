#' Import Ct values from QuantStudio & StepOne qPCR software
#'
#'import_ct() tidies the "Results" sheet of an .xls file produced by QuantStudio Design & Analysis Software v1.5.2 or StepOne Software v2.3.
#'
#' @param file character vector of length 1. Path to the .xls file containing the qPCR results to import.
#' @param from character vector of length 1. Indicates which qPCR software the input file originates from. Valid values are "QuantStudio" or "StepOne".
#'
#' @return A tibble with 4 columns: "well_position", "sample_id", "tar_nm", and "ct".
#'
#' @export
#'
#' @importFrom readxl read_excel
#' @importFrom stats setNames
#'
#' @examples
#'
import_ct <- function(file, from = "QuantStudio") {

  # Capture user's arguments

  file <- eval(expr = substitute(file), envir = parent.frame())

  from <- eval(expr = substitute(from), envir = parent.frame())

  # Stop execution in case of invalid input

  if (!requireNamespace("readxl", quietly = TRUE)) {
    warning("\nThe readxl package must be installed to use this function.\n")
    return(NULL)
  }

  if (!file.exists(file)) {
    stop("\nFile \"", file, "\" does not exist.\n")
  }

  if (length(from) != 1L || !from %in% c("QuantStudio", "StepOne")) {
    stop("\nfrom value must be \"QuantStudio\" or \"StepOne\".\n")
  }

  invalid_file <- unique(class(try(expr = match(x = "Well",
                                                table = suppressMessages(expr = read_excel(path = file,
                                                                                           sheet = "Results"))[[1]]),
                                   silent = T)) == "try-error")

  if (invalid_file) {
    stop("\nThe input file is not valid.\n")
  }

  skip <- match(x = "Well",
                table = suppressMessages(expr = read_excel(path = file,
                                                           sheet = "Results"))[[1]])

  missing_well_column <- ifelse(is.na(skip), T, F)

  if (missing_well_column) {
    stop("\nThe \"Results\" sheet in the input file does not contain the required \"Well\" column.\n")
  }

  # Read qPCR results file

  data <- suppressMessages(expr = read_excel(path = file,
                                             sheet = "Results",
                                             skip = skip))

  # Tidy QuantStudio input

  if (from == "QuantStudio") {

    if (!all(c("Well Position", "Sample Name", "Target Name", "CT") %in% colnames(data))) {
      if (all(c("Well", "Sample Name", "Target Name", paste0("C", "\u1D1B")) %in% colnames(data))) {
        stop("\nThe input file was generated by the StepOne software. Please provide from = \"StepOne\".\n")
      } else {
        stop("\nThe input file does not contain one or more required column(s).\n")
      }
    }

    data <- setNames(object = data[, c("Well Position", "Sample Name", "Target Name", "CT")],
                     nm = c("well_position", "sample_id", "tar_nm", "ct"))

  }

  # Tidy StepOne input

  if (from == "StepOne") {

    if (!all(c("Well", "Sample Name", "Target Name", paste0("C", "\u1D1B")) %in% colnames(data))) {
      if (all(c("Well Position", "Sample Name", "Target Name", "CT") %in% colnames(data))) {
        stop("\nThe input file was generated by the QuantStudio software. Please provide from = \"QuantStudio\".\n")
      } else {
        stop("\nThe input file does not contain one or more required column(s).\n")
      }
    }

    data <- setNames(object = data[c(1L:(nrow(data)-5L)), c("Well", "Sample Name", "Target Name", paste0("C", "\u1D1B"))],
                     nm = c("well_position", "sample_id", "tar_nm", "ct"))

  }

  # Tidy final output

  data$ct <- vapply(X = data$ct,
                    FUN = function(x) { ifelse(x == "Undetermined", NA_real_, as.double(x)) },
                    FUN.VALUE = vector(mode = "double", length = 1L))

  data <- data[order(data$sample_id, data$tar_nm), , drop = F]

  return(data)

}
