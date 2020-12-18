#
#     R-functions to calculate tests of goodness of fit of graduations
#
#     ============================================================
#
#     Chi-square test
#
#     Arguments: Z = residuals Z.x
#                Npar = number of fitted parameters
#
Chi.Square <- function(Z, Npar){
  Chis2 <- sum(Z^2)
  DF <- length(Z) - Npar
  Sig.Pr <- 1 - pchisq(Chis2, DF)
  return(list(Chis2 = Chis2, DF = DF, Sig.Pr = Sig.Pr))
}
#
#     ===========================================================
#
#     Standardized deviations test
#
#     Equal area test
#
#     Arguments: Z = residuals Z.x
#                N = Number of cells
#
Standard.Area <- function(Z, N){
  DF <- N-1
  Exp <- length(Z)/N
  CDF <- seq(0, 1, length = (N+1))
  Inner.Boundary <- qnorm(CDF[-c(1, length(CDF))])
  Boundary <- c(-20, Inner.Boundary, 20)
  Obs <- hist(Z, breaks = Boundary, plot = FALSE)$counts
  Chis2 <- sum( (Obs - Exp)^2/Exp)
  Sig.Pr <- 1 - pchisq(Chis2, DF)
  return(list(Boundary = round(Inner.Boundary, digits = 3), Obs = Obs, Exp = Exp,
              DF = DF, Chis2 = Chis2, Sig.Pr = Sig.Pr))
}
#
#     ==========================================================
#
#     Standardized deviations test
#
#     Equal width test
#
#     Arguments: Z = residuals Z.x
#                W = Width of internal cells
#
Standard.Width <- function(Z, W){
#
#   Calculate position of upper boundary and cell boundaries
#
  n <- length(Z)
  Upper <- qnorm(1 - 5/n)
  Upper <- max(seq(0, Upper, by = W))
  Inner.Boundary <- seq(-Upper, Upper, by = W)
  Boundary <- c(-20, Inner.Boundary, 20)
#
# Calculate DF and O_i
#
  DF <- length(Boundary) - 2
  Obs <- hist(Z, breaks = Boundary, plot = FALSE)$counts
#
# Calculate E_i
#
  CDF.p <- c(0, pnorm(Inner.Boundary), 1)  #  Cumulative probabilities
  CDF.f <- n * CDF.p                       #  Cumulative frequencies
  Exp <- diff(CDF.f)                       #  E_i
#
# Compute Chi^2 & sig prob
#
  Chis2 <- sum( (Obs - Exp)^2/Exp)
  Sig.Pr <- 1 - pchisq(Chis2, DF)
  return(list(Boundary = round(Inner.Boundary, digits = 3), Obs = Obs, Exp = Exp,
              DF = DF, Chis2 = Chis2, Sig.Pr = Sig.Pr))
}
#
#     ========================================================
#
#     Sign test
#
#     Argument: Z = residuals Z.x
#
Sign <- function(Z){
  n <- length(Z)
  Greater <- sum(Z >= 0)
  Less <- sum(Z < 0)
  if(Greater == Less) return("Number + is equal to number -")
  if(Greater > Less){
    Sig.Prob <- 2 * (1 - pbinom(Greater - 1, n, 0.5))
    return(list(N.plus = Greater, N.minus = Less, Sig.Prob = Sig.Prob))
  }
  if(Greater < Less){
    Sig.Prob <- 2 * pbinom(Greater, n, 0.5)
    return(list(N.plus = Greater, N.minus = Less, Sig.Prob = Sig.Prob))
  }
}
#
#     ========================================================
#
#     Change of sign test
#
#     Argument: Z = residuals Z.x
#
Change.Sign <- function(Z){
  Change <- 0
  N1 <- length(Z) - 1
  for(i in 1:N1) {
    if(Z[i] * Z[i+1] < 0) Change <- Change + 1
  }
  Sig.Pr <- pbinom(Change, N1, 0.5)
  return(list(N = N1+1, Change = Change, Sig.Pr = Sig.Pr))
}
#
#     =====================================================
#
#     Runs test
#
#     Arguments: Z = residuals Z.x
#
#     This function calls the subsiduary function Runs( )
#
Runs.test <- function(Z){
  Code <- Z; Code[ Z <= 0] = -1; Code[ Z > 0] = 1
  n1 <- sum(Code == 1); n2 <- sum(Code == -1)
  g <- 0; if(Code[1] > 0) g <- 1
  for(i in 1:(length(Code)-1)) if(Code[i] < Code[i+1]) g <- g + 1
  Sig.Pr <- Runs(n1, n2, g)
  return(list(n1 = n1, n2 = n2, g = g, Sig.Prob = Sig.Pr))
}
#
#     Subsiduary function to calculate significance probability
#
Runs <- function(n1, n2, g){
  Sig.Pr <- 0
  Denom <- choose(n1+n2, n1)
  for(i in 1:g){
    Num <- choose(n1-1,i-1)*choose(n2+1,i)
    Sig.Pr <- Sig.Pr + Num/Denom
  }
  Sig.Pr
}
#
#     ================================================
#
#     Runs test using a permutation test
#
#     Arguments: Z = residuals Z.x
#                n = number of samples - defaults to 1000 if n is not set
#     Calls the function Runs.test( )
#
Runs.Test.Perm <- function(Z, n = 1000){
  Code <- Z; Code[ Z <= 0] = -1; Code[ Z > 0] = 1
  g <- Runs.test(Z)$g
  Null.dist <- NULL
  for(i in 1:n) {
    Perm <- sample(Code)
    Runs <- 0; if(Perm[1] > 0) Runs <- 1
    for(i in 1:(length(Perm)-1)) if(Perm[i] < Perm[i+1]) Runs <- Runs + 1
    Null.dist <- c(Null.dist, Runs)
  }
  Sig.Pr <- sum(Null.dist <= g)/n
# Plot null distribution with observed g indicated
  Null.dist <- table(Null.dist)/n
  x <- as.numeric(names(Null.dist))
  par(mfrow = c(1,1))
  plot(x, Null.dist, type = "h", lwd = 2, col = "blue",
     main = "Null distribution of runs statistics",
     xlab = "Number of runs", ylab = "Probability")
  lines(g, Null.dist[paste(g)], type = "h", col = "red", lwd = 3)
  return(list(Runs = g, Null.dist = Null.dist, Sig.Pr = Sig.Pr))
}
#  
#     =================================================
#
#     Serial correlation test
#
#     Arguments: Z = residuals Z.x
#
Serial <- function(Z){
  Set1 <- Z[-length(Z)]; Set2 <- Z[-1]
  Corr <- cor(Set1, Set2)
  Sig.Pr <- 1 - pnorm(Corr, 0, sqrt(1/(length(Z)-1)))
  return(list(Serial = Corr, Sig.Pr = Sig.Pr))
}


