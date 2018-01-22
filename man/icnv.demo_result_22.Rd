\docType{data}
\name{icnv.demo_result_22}
\alias{icnv.demo_result_22}
\title{
    Example iCNV calling results.
}
\format{A matrix containing the coverage information for each nucleotide
    position of a ribosomal RNA 18S obtained in a structure probing experiment
    using DMS as a probe, where:
    \describe{
        \item{1st column}{control experimental replicate}
        \item{2nd column}{control experimental replicate}
        \item{3rd column}{control experimental replicate}
        \item{4th column}{treatment experimental replicate}
        \item{5th column}{treatment experimental replicate}
        \item{6th column}{treatment experimental replicate}
}}
\usage{
    covFile
}
\value{
    Coverage information per nucleotide.
}
\description{
    A matrix containing 1800x6 entries of per nucleotide coverage information
    from structure probing experiments.

    This matrix is provided as a reference for creating a
    \code{SummarizedExperiment} object for storing all structure probing data
    used by the \code{BUMHMM} package. See \code{\link[BUMHMM]{se}} for the
    example code.
}
\keyword{datasets}