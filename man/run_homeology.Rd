% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/predict.R
\name{run_homeology}
\alias{run_homeology}
\title{run_homeology}
\usage{
run_homeology(
  junctions,
  width,
  pad,
  thresh,
  stride,
  genome,
  cores,
  flip,
  bidirectional,
  annotate,
  savegMatrix
)
}
\arguments{
\item{junctions}{Path to the junctions file}

\item{width}{Width of the bins}

\item{pad}{Number of bases of padding around each sequence position (bin) to use when computing homeology, i.e. we then will be comparing 1 + 2*pad -mer sequences for edit distance}

\item{thresh}{String distance threshold for calling homeology in a bin}

\item{stride}{Distance in bases between consecutive bins in which we will be measuring homeology}

\item{genome}{Path to .2bit or ffTrack .rds containing genome sequence}

\item{cores}{How many cores to use}

\item{flip}{If flip = FALSE, homeology search for -/- and +/+ junctions is done between a sequence and its reverse complement}

\item{bidirectional}{Adding padding on both sides of each breakpoint (TRUE) or only in the direction of the fused side (FALSE)}

\item{annotate}{Annotate edges in gGraph object and save it in working directory}

\item{savegMatrix}{Save gMatrix object of edit distances}
}
\description{
function to run homeology
}
