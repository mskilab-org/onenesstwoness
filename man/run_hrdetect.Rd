% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/predict.R
\name{run_hrdetect}
\alias{run_hrdetect}
\title{run_hrdetect}
\usage{
run_hrdetect(
  snv,
  indel,
  jabba,
  sv,
  mask,
  hets,
  genome,
  ref,
  save = TRUE,
  outdir
)
}
\arguments{
\item{snv}{Path to the SNV file}

\item{indel}{Path to the indel file if not present in SNV}

\item{jabba}{Jabba rds file with $asegstats}

\item{sv}{Path to the SV file (vcf/rds)}

\item{mask}{Path to the mask file}

\item{hets}{Path to the hets file}

\item{genome}{enum of hg19/hg38/mm10/canFam3}

\item{ref}{Path to the reference genome}

\item{save}{Logical indicating whether to save the results}

\item{outdir}{Path to the output directory}
}
\description{
function to run HRDetect
}
