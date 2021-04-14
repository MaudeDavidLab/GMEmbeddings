\name{getSeqtab}
\alias{getSeqtab}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
getSeqtab
%%  ~~function to do ... ~~
}
\description{
Function reads in your seqtab file and puts it in a format that can be used as an argument in the EmbedAsvTable function.
}
\usage{
getSeqtab(seqtab)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{seqtab}{
seqtab should be a .csv file.

The seqtab argument should given as the file path and needs to be in quotes.
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
function (seqtab)
{
    seqtab <- utils::read.csv(seqtab)
    rownames(seqtab) <- seqtab$X
    rownames(seqtab) <- gsub("X", "", rownames(seqtab))
    seqtab <- seqtab[, 2:ncol(seqtab)]
    return(seqtab)
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
