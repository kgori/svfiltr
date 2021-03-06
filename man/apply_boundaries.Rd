% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/filters.R
\name{apply_boundaries}
\alias{apply_boundaries}
\title{Filter data frame so all samples selected by `index` have read support
greater than the boundary threshold, or 0}
\usage{
apply_boundaries(df, boundaries, index, lower_grace = 0, upper_grace = 0)
}
\arguments{
\item{df}{A dataframe}

\item{boundaries}{Estimated boundaries from estimated_boundaries}

\item{index}{Columns of df to examine}

\item{lower_grace}{A 'grace' number of reads: below this many reads an SV is considered absent}

\item{upper_grace}{A 'grace' number of reads: if an SV is this many reads below the boundary, it still counts as present}
}
\value{
Filtered dataframe
}
\description{
Filter data frame so all samples selected by `index` have read support
greater than the boundary threshold, or 0
}

