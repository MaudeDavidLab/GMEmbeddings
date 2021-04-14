\name{getEmbeddingMatrix}
\alias{getEmbeddingMatrix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
getEmbeddingMatrix
%%  ~~function to do ... ~~
}
\description{
Function reads in the embedding_file given by the user.
}
\usage{
getEmbeddingMatrix(embedding_file_name)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{embedding_file_name}{
This should given as the file path, and needs to be in quotes.

Ex: getEmbeddingMatrix("directory/embedding_file_name")
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
function (embedding_file_name)
{
    embedding_matrix <- read.csv(embedding_file_name, row.names = 1,
        sep = "")
    return(embedding_matrix)
  }

## To read in the example file that comes preinstalled with the package,
    type the following command:
    getEmbeddingMatrix(system.file("extdata", "embed_.07_100dim.txt", package = "GMEmbeddings"))

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
