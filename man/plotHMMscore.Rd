% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotHMMscore_function.R
\name{plotHMMscore}
\alias{plotHMMscore}
\title{Plot CNV inference score.}
\usage{
plotHMMscore(icnv_res, h = NULL, t = NULL, title = "score plot",
  output = NULL, col = "")
}
\arguments{
\item{icnv_res}{CNV inference result. Result from iCNV_detection() (i.e. iCNV_detection(...))}

\item{h}{start position of this plot. Default Start of the whole chromosome}

\item{t}{end position of this plot. Default End of the whole chromosome}

\item{title}{of this plot. Character value. Type character Default "score plot"}

\item{output}{generated from output_list_function. If it isn't null, only CNVs in output file will be highlighted. Default NULL}

\item{col}{Specify if would like to plot in DGV color scheme ('DGV',red for deletion, blue for duplication and grey for diploid) or default color scheme (blue for deletion, red for duplicatin and and green for diploid) Type character. Default ''}
}
\value{
void
}
\description{
Plot out CNV inference score. Each row is a sample, each column is a SNP or, 
exon (WES) or bin (WGS). Red color indicate score favor duplication whereas 
blue favor deletion.
}
\examples{
plotHMMscore(icnv_res0,h=21000000, t=22000000, title='my favorite subject')
plotHMMscore(icnv_res0,h=21000000, t=22000000, title='my favorite subject',col='DGV')
}
