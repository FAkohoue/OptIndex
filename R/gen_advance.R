#' Calculate Genetic Advance for a trait or combination of traits
#' @param phen_mat Phenotypic matrix value of desired characters
#' @param gen_mat Genotypic matrix value of desired characters
#' @param weight_mat Weight matrix value of desired characters
#' @param k Selection intensity based on the breeder's objective
#' @return Genetic advance of character or character combinations
#' @export
#'
#' @examples
#' weight_h <-read.table("Weight_varselect1.txt",h=TRUE, row.names=1)
#' pvar_cov <-read.table("pvar_varselect1.txt",h=TRUE, row.names=1)
#' gvar_cov <-read.table("gvar_varselect1.txt",h=TRUE, row.names=1)
#' (X <- as.matrix(Phenotyping))
#' (P <- as.matrix(pvar_cov))
#' (G <- as.matrix(gvar_cov))
#' (a1 <- as.matrix(weight_h))
#' gen.advance(phen_mat = P[1,1], gen_mat = G[1,1], weight_mat = weight[1,2], k = 1.271)
#'
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
