\name{muscle-package}
\alias{muscle-package}

\docType{package}
\title{
Multiple Sequence Alignment
}
\description{
MUSCLE performs multiple sequence alignments of nucleotide and amino acid sequences.
}
\details{
Details about the algorithm can be found on the MUSCLE website:\cr
\url{http://www.drive5.com/muscle/}
}

\author{
Algorithm: Robert C. Edgar

R port: Alex T. Kalinka \email{alex.t.kalinka@gmail.com}

}
\references{
Edgar, R.C. (2004) MUSCLE: multiple sequence alignment with high accuracy and high throughput. \emph{Nucleic Acids Res} \bold{32}, 1792-1797.
}

\seealso{
\code{\link{muscle}}, \code{\link{umax}}

}
\examples{
## Align sequences from a fasta file.
\dontrun{
aln <- muscle(seqs = "unaligned.fasta")
}

## Align sequences from a fasta object.
\dontrun{
aln <- muscle(seqs = umax)
}
}

