\name{EmbedAsvTable}
\alias{EmbedAsvTable}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
EmbedAsvTable

%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
EmbedAsvTable(seqtab, fasta_file, blast_hits, embedding_matrix)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{seqtab}{
The \code{seqtab} argument needs to be an object. If not already an object, use the \code{getSeqtab} function to change the \code{seqtab} file into something that the \code{EmbedAsvTable} function can use as an argument.}

  \item{fasta_file}{
The \code{fasta_file} argument needs to be a string. This string will contain the full file name of your fasta file.
}

  \item{blast_hits}{
  The \code{blast_hits} argument needs to be passed in as an object. If not an object, the \code{EmbedAsvTable} function will not work properly.
  }
  \item{embedding_matrix}{
  The \code{embedding_matrix} argument needs to be passed in as a string. This string should contain the full file name of your embedding matrix.
  }
}

\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

## The function is currently defined as
EmbedAsvTable <- function(seqtab, fasta_file_name, blast_hits, embedding_file_name){

  best_hits = getBestHits(blast_hits = blast_hits, id_thresh = 99)

  seqtab <- embedSeqtab(seqtab, fasta_file = fasta_file_name, best_hits = best_hits, embedding_file_name)

  seqtab
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
