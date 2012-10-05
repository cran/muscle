\name{read.fasta}
\alias{read.fasta}

\title{
Read a FASTA file into R
}
\description{
Reads sequence files in the FASTA format into an R object
}
\usage{
read.fasta(file)
}

\arguments{
  \item{file}{
A character string naming a fasta file.
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
\code{\link{write.fasta}}
}
\examples{
## Read in unaligned sequences.
\dontrun{
seqs <- read.fasta(file = "unaligned.fasta")
}
}

