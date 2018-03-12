\name{apt.cytoscan.process.batch}
\alias{apt.cytoscan.process.batch}
\title{Run apt.cytoscan.process() with multithreading}
\description{
  Same as apt.cytoscan.process() but in batch mode, with multithreading.
}
\usage{
  apt.cytoscan.process.batch(CEL.list.file = NULL, nthread = 1,
  cluster.type = "PSOCK", ...)
}
\arguments{
  \item{CEL.list.file}{A CEL filetab-separated file giving the file names of the CEL files and corresponding sample name (see Details section).}
  \item{nthread}{Number of threads for multithreading (see Details).}
  \item{cluster.type}{Type of multithreading cluster (see \code{?parallel::makeCluster()}).}
  \item{...}{Any parameter for the function \code{apt.cytoscan.process()}, except \code{CEL} and \code{samplename}.}
}
\details{
  \code{CEL.list.file} is a tab-separated text file containing 2 columns :
  \itemize{
    \item{cel_files : Name (and path) of the CEL file(s)}
    \item{SampleName : The output sample name}
  }
}
\author{
  Bastien Job
}

