#' @title Fits the beta distribution by method of moments.
#' @usage beta.mom(dat)
#' @param dat A vector contains the numbers that are assumed to follow a
#'           beta distribution.
#' @return Returns the estimation of alpha and beta.
#' @author Ning Leng
#' @examples beta.mom(rbeta(5,5,20))
beta.mom <- function(dat){
mu <- mean(dat)
sig <- var(dat)
ab <- mu*(1-mu)/sig -1
a <- ab*mu
b <- ab-a
out <- c(a,b)
}

