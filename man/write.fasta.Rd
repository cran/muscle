\name{write.fasta}
\alias{write.fasta}

\title{
Write a FASTA file
}
\description{
Writes a sequence object to a FASTA file format
}
\usage{
write.fasta(aln, file)
}

\arguments{
  \item{aln}{
An object of class \code{muscle} or \code{fasta}.
}
  \item{file}{
A character string naming a file to write to.
}

}

\value{
Writes a fasta file as a side-effect.
}

\author{
Alex T. Kalinka \email{alex.t.kalinka@gmail.com}
}


\seealso{
\code{\link{read.fasta}}
}
\examples{
## Write out aligned sequences.
\dontrun{
write.fasta(aln, file = "out.afa")
}
}

