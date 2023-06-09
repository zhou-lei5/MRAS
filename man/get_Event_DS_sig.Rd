% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_Event_DS_sig.R
\name{get_Event_DS_sig}
\alias{get_Event_DS_sig}
\title{This function provides a way to get significantly differentially spliced events from the differentially spliced event matrix.}
\usage{
get_Event_DS_sig(Events_DS, DS_pvalue = 0.05, DS_dPSI = 0.1)
}
\arguments{
\item{Events_DS}{The difference matrix of the clipping events, obtained by the MRAS function DS_matrix().}

\item{DS_pvalue}{Significance level of differential splicing analysis.}

\item{DS_dPSI}{Threshold for differential splicing.}
}
\value{
Events with significant differences clipped by threshold filtering in Events_DS.
}
\description{
This function provides a way to get significantly differentially spliced events from the differentially spliced event matrix.
}
