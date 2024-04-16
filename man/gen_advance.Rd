\name{gen_advance}
\alias{gen_advance}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Calculation of genetic advance for a trait or combination of traits
%%  ~~function to do ... ~~
}
\description{This function is designed to calculate the genetic advance for a specific trait or a combination of traits in a breeding program. Genetic advance is a measure of the expected improvement in the traits of the offspring, compared to the parents, as a result of selection based on certain criteria.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
gen.advance(phen_mat, gen_mat, weight_mat, k)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{phen_mat}{Phenotypic matrix value of desired characters
%%     ~~Describe \code{x} here~~
}
 \item{gen_mat}{Genotypic matrix value of desired characters
%%     ~~Describe \code{x} here~~
}
 \item{weight_mat}{Weight matrix value of desired characters
%%     ~~Describe \code{x} here~~
}
 \item{k}{Selection intensity based on the breeder's objective
%%     ~~Describe \code{x} here~~
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
\author{Félicien Akohoue
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{comb_indices_fa}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
#' weight_h <-read.table("Weight_varselect1.txt",h=TRUE, row.names=1)
#' pvar_cov <-read.table("pvar_varselect1.txt",h=TRUE, row.names=1)
#' gvar_cov <-read.table("gvar_varselect1.txt",h=TRUE, row.names=1)
#' (X <- as.matrix(Phenotyping))
#' (P <- as.matrix(pvar_cov))
#' (G <- as.matrix(gvar_cov))
#' (a1 <- as.matrix(weight_h))
#' gen.advance(phen_mat = P[1,1], gen_mat = G[1,1], weight_mat = weight[1,2], k = 1.271)
#'
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
gen_advance <- function(phen_mat, gen_mat, weight_mat, k){
  p <- as.matrix(phen_mat)
  g <- as.matrix(gen_mat)
  w <- as.matrix(weight_mat)
  cat("Calculating bmat...\n")
  bmat <- solve(phen_mat) %*% gen_mat %*% weight_mat
  cat("Calculating Genetic Advance...\n")
  GA <- k * t(bmat) %*% g %*% w / (t(bmat) %*% p %*% bmat)^0.5
  return(GA)
}
{
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
