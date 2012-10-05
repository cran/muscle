\name{muscle}
\alias{muscle}

\title{
Multiple Sequence Alignment
}
\description{
MUSCLE performs multiple sequence alignments of nucleotide or amino acid sequences.
}
\usage{
muscle(seqs, out = NULL, quiet = FALSE, ...)
}

\arguments{
  \item{seqs}{
Either a character string naming a FASTA file containing unaligned nucleotide or amino acid sequences, or an object of class \code{fasta}.
}
  \item{out}{
A character string naming the file to write aligned sequences to. If NULL then output is returned to the user as an R object of class \code{muscle}. Defaults to NULL.
}
  \item{quiet}{
Logical, whether to prevent printing of the progress of the algorithm to the screen. Defaults to FALSE.
}
  \item{\dots}{
Arguments (options and flags) for the MUSCLE algorithm (see Details).
}

}
\details{
Arguments for the MUSCLE algorithm are either options or flags. Options take a variety of values, but flags are always logicals (TRUE or FALSE).
All options and flags are explained in detail on the MUSCLE website: \url{http://www.drive5.com/muscle/muscle_userguide3.8.html}

}

\value{
Output is either written to a FASTA file (if \code{out} is not NULL), or saved to an object of class \code{muscle}, which is a list containing the following components:
\enumerate{
 \item{\code{seqs}  }{ A data frame of aligned sequences. The first column contains sequence names, and the second column contains aligned sequences.}
 \item{\code{length}  }{ The number of positions in the alignment,}
 \item{\code{num}  }{ The number of sequences.}
}
}
\references{
Edgar, R.C. (2004) MUSCLE: multiple sequence alignment with high accuracy and high throughput. \emph{Nucleic Acids Res} \bold{32}, 1792-1797.
}
\author{
Algorithm by Robert C. Edgar. Ported into R by Alex T. Kalinka \email{alex.t.kalinka@gmail.com}
}
\note{
For further details see the MUSCLE website:
\url{http://www.drive5.com/muscle/}
}


\seealso{
\code{\link{print.muscle}}, \code{\link{read.fasta}}, \code{\link{write.fasta}}, \code{\link{umax}}
}
\examples{
## Align sequences from a fasta file.
\dontrun{
aln <- muscle(seqs = "unaligned.fasta")
}

## Align sequences from a fasta object.
aln <- muscle(seqs = umax)

## Align sequences from a fasta object and write them to a fasta file.
\dontrun{
aln <- muscle(seqs = umax, out = "out.afa")
}

## Switch on the diags flag.
aln <- muscle(seqs = umax, diags = TRUE)

}

