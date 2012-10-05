\name{print.muscle}
\alias{print.muscle}

\title{
Print Alignment Regions
}
\description{
This function prints aligned sequences for specific locations and sequences.
}
\usage{
\method{print}{muscle}(x, ...)
}

\arguments{
  \item{x}{
An object of class \code{muscle}.
}
  \item{\dots}{
Alignment positions and/or sequences (see Details).
}
}

\details{
Additional arguments for print:
\enumerate{
 \item \code{from}: the first position (integer)
 \item \code{to}: the last position (integer)
 \item \code{num}: the number of sequences to display (integer)
 \item \code{seqs}: specific sequences to display (character vector).
}
}

\value{
Prints aligned sequences to the screen.
}

\author{
Alex T. Kalinka \email{alex.t.kalinka@gmail.com}
}

\examples{
## Align sequences.
aln <- muscle(seqs = umax)

## Print the beginning of the alignment to the screen.
print(aln)

## Print a longer stretch of the alignment.
print(aln, from = 1, to = 179)

## Print the last 150 positions in the alignment.
num <- aln$length
print(aln, from = num-150, to = num)

## Print a portion of the alignment for the first 5 sequences.
print(aln, from = 156, to = 289, num = 5)

## Print a portion of the alignment for a specific set of sequences.
print(aln, from = 156, to = 289, seqs = c("Sus_scrofa", "Homo_sapiens", "Felis_catus"))

## Print interleaved aligned sequences to a file.
\dontrun{
sink(file = "aln.il", type = "output")
print(aln, from = 1, to = aln$length)
sink()
}

}

