\name{comb_indices_fa}
\alias{comb_indices_fa}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Construction of selection indices based on a number of character groupings
%%  ~~function to do ... ~~
}
\description{construct selection indices based on a set of input parameters related to breeding values, variance-covariance matrices, weights, genetic advances, desired gains, heritability, and selection intensity.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
comb.indices_fa(BV, ncomb = 1, matp, matg, wmat, wcol = 1, GA, wg, d, h2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{BV}{Matrix of breeding values of genotypes for all traits
%%     ~~Describe \code{x} here~~
}
  \item{ncomb}{Number of Characters/Traits group
%%     ~~Describe \code{x} here~~
}
  \item{matp}{Phenotypic Variance-Covariance Matrix
%%     ~~Describe \code{x} here~~
}
  \item{matg}{Genotypic Variance-Covariance Matrix
%%     ~~Describe \code{x} here~~
}
  \item{wmat}{Weight Matrix
%%     ~~Describe \code{x} here~~
}
  \item{wcol}{Weight column number in case more than one weights, by default its 1
%%     ~~Describe \code{x} here~~
}
  \item{GA}{Genetic Advance of comparative Character/Trait i.e. Yield (Optional argument)
%%     ~~Describe \code{x} here~~
}
  \item{wg}{Subset of Genotypic Variance-Covariance Matrix
%%     ~~Describe \code{x} here~~
}
  \item{d}{Vector of desired gains of traits
%%     ~~Describe \code{x} here~~
}
  \item{h2}{Heritability of the trait of interest
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

\seealso{\code{\link{gen_advance}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
#' NG <- 408
#' NT <- 7
#' k <- 1.271
#' tr <- 1
#' weight_h <-read.table("Weight_varselect1.txt",h=TRUE, row.names=1)
#' pvar_cov <-read.table("pvar_varselect1.txt",h=TRUE, row.names=1)
#' gvar_cov <-read.table("gvar_varselect1.txt",h=TRUE, row.names=1)
#' (X <- as.matrix(Phenotyping))
#' (P <- as.matrix(pvar_cov))
#' (G <- as.matrix(gvar_cov))
#' (a1 <- as.matrix(weight_h))
#' (a1mat <- subset(a1, select = -Equal))
#' (a1mat <- as.matrix(a1mat))
#' (wg <- G[tr, tr])
#' H <- subset(weight_h, select = -Equal)
#' (h2 <- H[tr, 1])
#' g1 <- G[, tr]
#' g <- data.matrix(g1)
#' d <- sqrt(diag(G))
#'
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
comb_indices_fa <- function(BV, ncomb, matp, matg, wmat, wcol = 1, GA, wg, d, h2, k){
  cat("Starting the index calculation...\n")

    selection.index <- function(ID, matp, matg, wmat, GA){
    cat("Calculating selection index for ID:", ID, "\n")
    ID = toString(ID)
    p <- as.matrix(matp)
    g <- as.matrix(matg)
    w <- as.matrix(wmat)
    bmatI <- solve(p) %*% g %*% w
    G <- k * t(bmatI) %*% g %*% w / (t(bmatI) %*% p %*% bmatI)^0.5
    RE <- if(missing(GA)) (G/G) else (G/GA)
    result <- list("ID" = ID, "b" = matrix(round(bmatI,4), nrow = 1), "GA" = round(G,4), "RE" = round(RE,4))
    return(data.frame(result))
  }

    ncolmn <- ncol(matp)
    comb <- t(combn(ncolmn, ncomb))
    indices <- list()
    for (i in 1:nrow(comb)) {
    ID <- paste0(comb[i,], collapse = "")
    cat("Processing combination:", ID, "\n")
    indices[[i]] <- selection.index(ID,
                                    matp = matp[comb[i,], comb[i,]],
                                    matg = matg[comb[i,], comb[i,]],
                                    wmat = wmat[comb[i,], wcol], GA = GA)
  }

    df <- do.call(rbind.data.frame, indices)
    df$Rank <- rank(-df$RE, ties.method = "max")
    df$ID <- rownames(matp)
    rownames(df) <- df$ID
    df <- subset(df, select = -ID)

    cat("Finished calculating selection indices for all combinations.\n")

    bmat <- subset(df, select = b)

    if (!is.matrix(X)) {
      X <- as.matrix(X)
  }

    if (!is.matrix(bmat)) {
    bmat <- as.matrix(bmat)
  }

    # Check for NA values
    if (any(is.na(X))) {
    stop("NA values found in X")
  }

   if (any(is.na(bmat))) {
    stop("NA values found in bmat")
  }

    cat("Dimensions of X:", dim(X), "\n")
    cat("Dimensions of bmat:", dim(bmat), "\n")

    # Optimum index (Smith 1936)
    I_O <- X%*%bmat
    cat("Calculated Optimum index following Smith 1936\n")
    RHI_O <- (sqrt(t(bmat)%*%matp%*%bmat))/(sqrt(t(wmat)%*%G%*%wmat))
    GA_total_O <- k%*%RHI_O%*%(sqrt(t(wmat)%*%matg%*%wmat))
    GA_O <- subset(df, select = GA)
    colnames(GA_O) <- c("GA_O")
    rG_O <- (t(bmat)%*%g)/(sqrt(wg%*%(t(bmat)%*%matp%*%bmat)))
    RE_O <- rG_O/h2
    CV_O <- sqrt(t(bmat) %*% matp %*% bmat) / mean(bmat) * 100

    # Base index (Brim et al. 1959)
    b_B<-a1mat
    I_B <- X%*%b_B
    cat("Calculated base index following Brim et al. 1959\n")
    RHI_B <- (sqrt(t(a1mat)%*%matg%*%a1mat))/(sqrt(t(a1mat)%*%matp%*%a1mat))
    deltaH_B <- k%*%sqrt(t(a1mat)%*%matp%*%a1mat)
    delta_B <- (k*matg%*%a1mat)%*%solve(sqrt(t(a1mat)%*%matp%*%a1mat))
    rG_B <- (t(a1mat)%*%g)/(sqrt(wg%*%(t(a1mat)%*%matp%*%a1mat)))
    RE_B <- rG_B/h2
    CV_B <- sqrt(t(a1mat)%*%matp%*%a1mat)/mean(I_B)*100

    #Pesek and Baker index (Pesek and Baker  1969)
    b_PB <- solve(matg)%*%d
    I_PB <- X%*%b_PB
    cat("Calculated Selection index following Pesek and Baker 1969\n")
    RHI_PB <- (sqrt(t(b_PB)%*%matp%*%b_PB))/(sqrt(t(d)%*%matg%*%d))
    deltaH_PB <- k%*%RHI_PB%*%(sqrt(t(d)%*%matg%*%d))
    delta_PB <- (k*matg%*%b_PB)%*%solve(sqrt(t(b_PB)%*%matp%*%b_PB))
    rG_PB <- (t(b_PB)%*%g)/(sqrt(wg%*%(t(b_PB)%*%matp%*%b_PB)))
    RE_PB <- rG_PB/h2
    CV_PB <- sqrt (t(b_PB)%*%P%*%b_PB)/mean(I_PB)*100

    #Rank of indices
    rO <- rank (-I_O, ties.method = "max")
    mrank_O <- cbind(I_O,rO)
    colnames(mrank_O) <- c("Optimum_Index_Smith", "Gen_rank_Smith")
    rB <- rank (-I_B, ties.method = "max")
    mrank_B <- cbind(I_B,rB)
    colnames(mrank_B) <- c("Index_Base", "Gen_rank_Base")
    rPB <- rank (-I_PB, ties.method = "max")
    mrank_PB <- cbind(I_PB,rPB)
    colnames(mrank_PB) <- c("Index_Pesek_and_Baker", "Gen_rank_Pesek_and_Baker")

    out1 <- cbind(bmat,GA_O,b_B,delta_B,b_PB,delta_PB)
    colnames(out1) <- c("b_O", "GA_O", "b_B", "GA_B","b_PB", "GA_PB")

    # Correlations between ranks
    out2 <- cor(cbind(rO,rB,rPB))

    out3 <- cbind(RHI_O,GA_total_O,RE_O,CV_O, RHI_B,deltaH_B,RE_B,CV_B, RHI_PB,deltaH_PB,RE_PB,CV_PB)
    colnames(out3) <- c("RHI_O", "GA_total_O", "RE_O", "CV_O", "RHI_B", "deltaH_B", "RE_B", "CV_B", "RHI_PB", "deltaH_PB", "RE_PB", "CV_PB")


    out4 <- cbind(mrank_O,mrank_B,mrank_PB)

    # Create the Excel file
    wb <- createWorkbook()
    # Add each data frame to a separate sheet in the Excel file
    addWorksheet(wb, "Index_coefs")
    writeData(wb, sheet = "Index_coefs", x = out1, rowNames = TRUE)
    addWorksheet(wb, "Rank_Correlations")
    writeData(wb, sheet = "Rank_Correlations", x = out2, rowNames = TRUE)
    addWorksheet(wb, "Index_Evaluation")
    writeData(wb, sheet = "Index_Evaluation", x = out3, rowNames = TRUE)
    addWorksheet(wb, "Index_and_Ranking")
    writeData(wb, sheet = "Index_and_Ranking", x = out4, rowNames = TRUE)
    # Save the Excel file
    saveWorkbook(wb, file = "Selection_Indices_varselect.xlsx", overwrite = TRUE)

    cat("Excel file created and saved successfully\n")

    return(df)
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
