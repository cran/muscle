\name{XStringSet2fasta}
\alias{XStringSet2fasta}

\title{
Convert XStringSet objects into fasta objects
}
\description{
Converts an XStringSet object (DNAStringSet, RNAStringSet, AAStringSet) into a fasta object.
}
\usage{
XStringSet2fasta(stringset, names)
}

\arguments{
  \item{stringset}{
An XStringSet object (DNAStringSet, RNAStringSet, AAStringSet) from the Biostrings package (Bioconductor).
}

  \item{names}{
A character string vector giving names for each of the sequences in the XStringSet object.
}

}

\value{
An object of class \code{fasta}, which is a list containing the following elements:
\enumerate{
 \item{\code{seqs} }{A data frame in which the first column contains sequence names and the second column contains sequences.}
 \item{\code{num} }{The number of sequences.}
}
}

\author{
Alex T. Kalinka \email{alex.t.kalinka@gmail.com}
}

\seealso{
\code{\link{write.fasta}}, \code{\link{read.fasta}}
}
\examples{
## Convert a DNAStringSet object to a fasta object.
\dontrun{
xx <- DNAStringSet(c("AGCTTAGCT","AGTTAGC","CGTTAGTTA"))
ff <- XStringSet2fasta(xx, names = 1:3)
}
}

