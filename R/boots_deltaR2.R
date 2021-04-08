##' Bootstrapping procedure to obtain p-value for differential R2.
##'
##' Bootstrapping procedure to obtain p-value for differential R2.
##' @title Likelihood method to obtain p-value for differential R2.
##' @param tt1 Time vector of condition 1
##' @param yy1 Expression vector of condition 1
##' @param tt2 Time vector of condition 2
##' @param yy2 Expression vector of condition 2
##' @param period Period of the since curve. Default is 24.
##' @param B number of Bootstrapping
##' @param seed random seed to reproduce the result
##' @param method Bootstrapping method, either normal or percentile.
##' @return Bootstrapping p-value. 
##' Formula 1: \eqn{yy = amp \times sin(2\pi/period \times (phase + tt)) + offset}
##' Formula 2: \eqn{yy = A \times sin(2\pi/period \times tt) + B * cos(2*pi/period * tt) + offset}
##' @author Caleb
##' @export
##' @examples
##' set.seed(32608)
##' n <- 50
##' tt1 <- runif(n,0,24) 
##' Amp1 <- 2
##' Phase1 <- 6
##' Offset1 <- 3
##' yy1 <- Amp1 * sin(2*pi/24 * (tt1 + Phase1)) + Offset1 + rnorm(n,0,1)
##' tt2 <- runif(n,0,24) 
##' Amp2 <- 3
##' Phase2 <- 5
##' Offset2 <- 2
##' yy2 <- Amp2 * sin(2*pi/24 * (tt2 + Phase2)) + Offset2 + rnorm(n,0,1)
##' boots_deltaR2(tt1, yy1, tt2, yy2)

boots_deltaR2 <- function(tt1, yy1, tt2, yy2, period = 24, B=1000, seed = 15213, method = "normal"){
	
	n1 <- length(tt1)
	stopifnot(n1 == length(yy1))
	n2 <- length(tt2)
	stopifnot(length(tt2) == length(yy2))
	
	#period <- 24
	w <- 2*pi/period
	
	fit1 <- fitSinCurve(tt1, yy1, period = 24)
	fit2 <- fitSinCurve(tt2, yy2, period = 24)
	deltaR2_obs <- fit2$R2 - fit1$R2
	
	set.seed(seed)
	deltaR_boots <- numeric(B)
	
	for(b in 1:B){
		bindex1 <- sample(n1, replace = TRUE)
		bindex2 <- sample(n2, replace = TRUE)
		
		bfit1 <- fitSinCurve(tt1[bindex1], yy1[bindex1], period = 24)
		bfit2 <- fitSinCurve(tt2[bindex2], yy2[bindex2], period = 24)
		
		deltaR_boots[b] <- bfit2$R2 - bfit1$R2
	}
	
	if(method == "normal"){
		Z <- deltaR2_obs / sd(deltaR_boots)
		pvalue <- min(pnorm(abs(Z), lower.tail = FALSE) * 2, 1)
	} else if (method == "percentile"){
		amean <- mean(deltaR_boots > 0)
		bmean <- mean(deltaR_boots < 0)
		pvalue <- min(amean, bmean) * 2
	} else {
		stop("no such method")
	}	
	pvalue
}



