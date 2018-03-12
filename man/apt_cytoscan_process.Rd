\name{apt.cytoscan.process}
\alias{apt.cytoscan.process}
\title{Normalize CytoScanHD_Array arrays using Affymetrix Power Tools.}
\description{
  This function performs the normalization of CEL file list for Affymetrix CytoScanHD_Array microarrays using the binary apt-copynumber-cyto-ssa (v2.4.0 2015-08-28) from Affymetrix Power Tools.
}
\usage{
apt.cytoscan.process(CEL = NULL, samplename = NULL, dual.norm = FALSE,
normal.diploid = FALSE, out.dir = getwd(), temp.files.keep = FALSE,
force.OS = NULL, apt.build = "na33.r4")
}
\arguments{
\item{CEL}{A CEL file.}
\item{samplename}{Name of the analyzed sample.}
\item{dual.norm}{Perform a dual-channel normalization (recommended for a better wave-effect bias normalization).}
\item{normal.diploid}{Consider the sample as normal-diploid (not recommended for tumors).}
\item{out.dir}{The path to output results.}
\item{temp.files.keep}{Keep the (numerous) intermediate files.}
\item{force.OS}{Specifiy the OS type of the execution machine (see Details).}
\item{apt.build}{Affymetrix genome build version.}
}
\details{
  \code{CEL.list.file} is a tab-separated text file containing 2 columns :
  \itemize{
    \item{cel_files : Name (and path) of the CEL file(s)}
    \item{SampleName : The output sample name}
  }
  Output files (all begin with the sample name as a prefix) :
  \itemize{
    \item{'*.OSCHP' : Contains the processed data. This format actually internaly is HDF5 [TO UPDATE!!].}
  }
  \code{force.OS} is \code{NULL} by default, but can take any of these values :
  \itemize{
    \item{linux}
    \item{windows}
    \item{osx}
  }
}
\author{
Bastien Job
}

