% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/import_ct.R
\name{import_ct}
\alias{import_ct}
\title{Import Ct values from QuantStudio & StepOne qPCR software}
\usage{
import_ct(file, from = "QuantStudio")
}
\arguments{
\item{file}{character vector of length 1. Path to the .xls file containing the qPCR results to import.}

\item{from}{character vector of length 1. Indicates which qPCR software the input file originates from. Valid values are "QuantStudio" or "StepOne".}
}
\value{
A tibble with 4 columns: "well_position", "sample_id", "tar_nm", and "ct".
}
\description{
import_ct() tidies the "Results" sheet of an .xls file produced by QuantStudio Design & Analysis Software v1.5.2 or StepOne Software v2.3.
}
