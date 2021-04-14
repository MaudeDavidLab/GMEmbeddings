\name{getFasta}
\alias{getFasta}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
getFasta
%%  ~~function to do ... ~~
}
\description{
Function reads in your fasta file puts it in a format that can be used as an argument in the EmbedAsvTable function.
}
\usage{
getFasta(fasta_file_name)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fasta_file_name}{
The fasta_file_name should be given as the file path to your .fasta file and should be in quotes.
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
Christine Tataru
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
function (fasta_file_name)
{
    fasta_file <- getExampleFasta(fasta_file_name)
    fasta_df <- getFastaDF(fasta_file)
    return(fasta_df)
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
