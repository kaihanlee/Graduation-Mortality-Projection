#
#     R-code to generate data for Assignment.  This code MUST be
#     placed at the start of your own R-script.  You must edit
#     the argument to set.seed( ) function to fit your own
#     registration number
#
set.seed(3031)    ##### EDIT TO THE LAST 4 DIGITS OF YOUR STUDENT ID ######
#
#     Read in CMI data and load R-funcions
#
#     NB: The files CMI_read.r, the CMI data files and R-function file
#         Test_GoF.r must be in the appropriate directory
#
source("CMI_read.r")
source("Test_GoF.r")
Year_choice <- Year[26:36]
My_Year <- sample(Year_choice, 1); My_Year  ## Year for assigmnent
AGE <- Age[Age > 39]

## Data for Part 1

DTH_MY <- Dth[Age > 39, Year == My_Year]
EXP_MY <- Exp[Age > 39, Year == My_Year]

## Data for Part 2

DTH <- Dth[ Age>=40 & Age <= 90 , Year>=1950 & Year <=My_Year ]
EXP <- Exp[ Age>=40 & Age <= 90 , Year>=1950 & Year <=My_Year ]

######### Solutions below here  ###########
###1.1####
mu <- DTH_MY/EXP_MY   #maximum likelihood esimate of force of mortality
obs <- log(mu)    #log of mortality

#generalized linear model with poisson errors and log link
Gomp.glm <- glm(DTH_MY ~ offset(log(EXP_MY)) + AGE, family=poisson);summary(Gomp.glm)
Gomp.glm$coef     #coefficients of Gompertz GLM

#plot Gompertz generalized linear model with poisson errors and log link
plot(obs ~ AGE, xlab = "Age", ylab = "log(mortality)", main = "Gompertz law: CMI data, ages 40 to 80, year 1983")
lines(AGE, Gomp.glm$lin - log(EXP_MY), col="red") #Add Gompertz line using the linear predictor
lines(AGE, Gomp.glm$coef[1] + Gomp.glm$coef[2] * AGE, col="red") #Add Gompertz line using fitted coefficients
lines(AGE, log(Gomp.glm$fit / EXP_MY), col="red")  #Add Gompertz line using fitted values
#we can see that all three Gompertz lines have the same equation
#equation of Gompertz line
text(55, -5.5, sprintf("log(mu) =  %.3f + %.3f * AGE", Gomp.glm$coef[1],Gomp.glm$coef[2]), adj = 0, cex = 1.2)
text(40, -3.2, "Poisson GLM", adj = 0, cex = 1.2) #Poisson GLM

###1.2####
z.i <- resid(Gomp.glm, type = "pear")   #pearson residuals

#plot pearson residuals against age
plot(z.i ~ AGE, ylab = "residuals", 
     xlab = "Age",main = "Residuals in Gompertz model (Poisson errors)")
points(AGE, resid(Gomp.glm), pch = 2,col="red") #points of deviance residuals
abline(h=c(-1.96,0,1.96), col="blue",lty=c(2,1,2))  #plot line at h=0 and 95% confidence intervals
legend(83, 6, legend = c("Pearson", "Deviance"), pch = c(1,2),col=c("black","red")) #legend

anova(Gomp.glm,test="Chisq")  #analysis of variance for Gompertz model

###Statistical tests####
#Chi-square test
# Arguments: Z = residuals z.i
#            Npar = number of fitted parameters
Chi.Square <- function(Z, Npar){
  Chis2 <- sum(Z^2)
  DF <- length(Z) - Npar
  Sig.Pr <- 1 - pchisq(Chis2, DF)
  return(list(Chis2 = Chis2, DF = DF, Sig.Pr = Sig.Pr))
}
Chi.Square(z.i,2)

#Standardized deviations test (Equal Area)
# Arguments: Z = residuals z.i
#            Npar = number of cells
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
Standard.Area(z.i,10)

#Standardized deviations test (Equal Width)
# Arguments: Z = residuals z.i
#            Npar = width of internal cells
Standard.Width <- function(Z, W){
  
  #Calculate position of upper boundary and cell boundaries
  n <- length(Z)
  Upper <- qnorm(1 - 5/n)
  Upper <- max(seq(0, Upper, by = W))
  Inner.Boundary <- seq(-Upper, Upper, by = W)
  Boundary <- c(-20, Inner.Boundary, 20)
  
  #Calculate DF and O_i
  DF <- length(Boundary) - 2
  Obs <- hist(Z, breaks = Boundary, plot = FALSE)$counts
  
  #Calculate E_i
  CDF.p <- c(0, pnorm(Inner.Boundary), 1)  #Cumulative probabilities
  CDF.f <- n * CDF.p                       #Cumulative frequencies
  Exp <- diff(CDF.f)                       #E_i

  #Compute Chi^2 & sig prob
  Chis2 <- sum( (Obs - Exp)^2/Exp)
  Sig.Pr <- 1 - pchisq(Chis2, DF)
  return(list(Boundary = round(Inner.Boundary, digits = 3), Obs = Obs, Exp = Exp,
              DF = DF, Chis2 = Chis2, Sig.Pr = Sig.Pr))
}
Standard.Width(z.i,0.5)

#Sign test
# Arguments: Z = residuals z.i
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
Sign(z.i)

#Change of sign test
# Arguments: Z = residuals z.i
Change.Sign <- function(Z){
  Change <- 0
  N1 <- length(Z) - 1
  for(i in 1:N1) {
    if(Z[i] * Z[i+1] < 0) Change <- Change + 1
  }
  Sig.Pr <- pbinom(Change, N1, 0.5)
  return(list(N = N1+1, Change = Change, Sig.Pr = Sig.Pr))
}
Change.Sign(z.i)

#Subsiduary function to calculate significance probability
Runs <- function(n1, n2, g){
  Sig.Pr <- 0
  Denom <- choose(n1+n2, n1)
  for(i in 1:g){
    Num <- choose(n1-1,i-1)*choose(n2+1,i)
    Sig.Pr <- Sig.Pr + Num/Denom
  }
  Sig.Pr
}

#Runs test
# Arguments: Z = residuals z.i
Runs.test <- function(Z){
  Code <- Z; Code[ Z <= 0] = -1; Code[ Z > 0] = 1
  n1 <- sum(Code == 1); n2 <- sum(Code == -1)
  g <- 0; if(Code[1] > 0) g <- 1
  for(i in 1:(length(Code)-1)) if(Code[i] < Code[i+1]) g <- g + 1
  Sig.Pr <- Runs(n1, n2, g)
  return(list(n1 = n1, n2 = n2, g = g, Sig.Prob = Sig.Pr))
}
Runs.test(z.i)

#Runs test using a permutation test
# Arguments: Z = residuals Z.x
#            n = number of samples - defaults to 1000 if n is not set
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
Runs.Test.Perm(z.i)

#Serial correlation test
# Arguments: Z = residuals z.i
#            Npar = number of cells
library(mgcv)
Serial <- function(Z){
  Set1 <- Z[-length(Z)]; Set2 <- Z[-1]
  Corr <- cor(Set1, Set2)
  Sig.Pr <- 1 - pnorm(Corr, 0, sqrt(1/(length(Z)-1)))
  plot(Z ~ AGE, ylab = "Residual", main = "Pearson Residual plot for Poisson model for mu")
  lines(AGE, gam(Z ~ s(AGE))$fitted, col = "blue", lwd = 2)
  return(list(Serial = Corr, Sig.Pr = Sig.Pr))
}
Serial(z.i)

sum(DTH_MY-Gomp.glm$fit)/sqrt(sum(Gomp.glm$fit))  #cumulative deviation

###2.1####
#Load gnm library
library("gnm")

MU <- DTH/EXP     #maximum likelihood estimate of force of mortality
Obs <- log(MU)    #log of mortality

#Convert to vectors
Dth.V <- c(DTH)
Exp.V <- c(EXP)

#Age and Year as factors
AGE <- 40:90
YEAR <- 1950:My_Year
Age.F <- factor(rep(AGE , ncol(DTH)))
Year.F <- factor(rep(YEAR , each=nrow(DTH)))  

#Fit model
LC.Model <- gnm( Dth.V ~ -1 + Age.F + Mult(Age.F,Year.F) , offset=log(Exp.V) , family=poisson)
Alpha.gnm <- LC.Model$coefficients[1:length(AGE)] 
Beta.gnm <- LC.Model$coefficients[(length(AGE)+1):(2*length(AGE))] 
Kappa.gnm <- LC.Model$coefficients[(2*(length(AGE))+1):((2*(length(AGE))+1)+(My_Year-1950))] 

#Satisfy identifiability constraints
Kappa.m <- mean(Kappa.gnm) 
Beta.m <- mean(Beta.gnm) 
Alpha.hat <- Alpha.gnm + Kappa.m * Beta.gnm 
Beta.hat <- Beta.gnm / (nrow(DTH) * Beta.m) 
Kappa.hat <- nrow(DTH) * Beta.m * (Kappa.gnm - Kappa.m)

#Fitted values
Fitted.M.hat = Alpha.hat + Beta.hat %*% t(Kappa.hat)

#function to plot Lee-Carter model
PlotLeeCarter <- function( AGE , YEAR , Obs , Alpha.hat , Beta.hat , Kappa.hat , Fitted.M.hat )
{
  par(mfrow = c(2,2), mar = c(4.5,4.5,1,1))
  #plot for alpha against age
  plot(AGE, Alpha.hat, xlab = "Age", ylab = expression(alpha), cex = 0.5, pch = 16)
  #plot for beta against age
  plot(AGE, Beta.hat, xlab = "Age", ylab = expression(beta), cex = 0.5, pch = 16)
  #plot for kappa against age
  plot(YEAR, Kappa.hat, xlab = "Year", ylab = expression(kappa), cex = 0.5, pch = 16)
    
  Age.Plot = 65  # Set plotting age
  Row = Age.Plot - min(AGE) + 1
  #plot Lee-Carter model
  plot(YEAR + 0.5, Obs[Row, ], xlab = "Year", ylab =
         expression(paste("log(",MU,")")), cex = 0.5, pch = 16)
  #plot best fit line
  lines(YEAR + 0.5, Fitted.M.hat[Row, ], type="l")
  legend("topright", legend = c("Data", "LC"), pch = c(16, -1), lty = c(-1, 1),
         bty = "n")     #legend
  legend("bottomleft", legend = paste("Age", Age.Plot), bty = "n")  #legend
}

#plot results
PlotLeeCarter( AGE , YEAR , Obs , Alpha.hat , Beta.hat , Kappa.hat , Fitted.M.hat )

1-pchisq(3062.079,1600) #chi-squared test for Lee-Carter model

#forecast kappa with drift model
library("astsa")  #nstall packages("astsa")
N.Ahead = 20    #20 years projection
par(mfrow=c(1,1))   #fit one graph in one interface
Kappa.for = sarima.for(Kappa.hat, n.ahead = N.Ahead, p=0, d=1, q=0) #kappa drift model
n.y <- ncol(DTH)    #count years

#Prediction error
Central = Kappa.for$pred    #forecasted kappa
SE.Pred = Kappa.for$se      #standard error of kappa
Z = qnorm(0.975)      #95% CI (two-tailed)
Kappa.Up.Pred = Central + Z*SE.Pred     #Upper CI
Kappa.Dn.Pred = Central - Z*SE.Pred     #Lower CI
Forecast = Alpha.hat + Beta.hat %*% t(Central)    #fit into Lee-Carter model
Forecast.Up.Pred = Alpha.hat + Beta.hat %*% t(Kappa.Up.Pred)  #Upper CI
Forecast.Dn.Pred = Alpha.hat + Beta.hat %*% t(Kappa.Dn.Pred)  #Lower CI
Range = (My_Year+1):(My_Year+20)  #forecast range

#Graphical output at age 65
LCForePred <- function(){   #Lee-Carter forecast model
  Plot.Age = 65
  Plot.Row = Plot.Age - min(AGE) + 1
  par(mfrow = c(1,1))
  par(mar=c(4.2, 4, 1, 0.5), mgp=c(3, 1, 0), las=1, cex=1)
  #plot Lee-Carter model for existing points
  plot(YEAR + 0.5, Obs[Plot.Row, ], axes = FALSE, xlab = "Year", ylab = "log(mortality)",
       xlim = c(1950,2003), ylim = c(min(Forecast.Dn.Pred[Plot.Row, ]), max(Obs[Plot.Row, ]) ), 
       main="Lee-Carter model fitted to males, CMI, ages 40-90, from 1950 to 1983")
  axis(1,las = 1, at = seq(1950, 2003, by = 10), tcl = -0.4)
  axis(2, seq(-6,-2, by = 0.5), tcl = -0.4)
  #add forecasted log of mortality and its confidence interval
  lines(Range + 0.5, Forecast[Plot.Row, ], lwd = 2)
  lines(Range + 0.5, Forecast.Up.Pred[Plot.Row, ], lty = 2, lwd = 2)
  lines(Range + 0.5, Forecast.Dn.Pred[Plot.Row, ], lty = 2, lwd = 2)
}
LCForePred()  #plot graph
legend("bottomleft", legend = c("Observed", "Central forecast", "95% CI Prediction error"),
       lty = c(-1, 1, 2, 2), pch = c(1,-1,-1,-1), lwd = 2, bty = "n")   #legend

#Add actual mortality 1984-2003
Dth.obs <- Dth[ Age==65 , Year>=My_Year+1 & Year <=My_Year+20 ]   #obtain actual deaths
Exp.obs <- Exp[ Age==65 , Year>=My_Year+1 & Year <=My_Year+20 ]   #obtain actual exposed
Obs.obs <- log( Dth.obs/Exp.obs )   #log of actual mortality

LCForePred()  #plot graph
points( (My_Year+1):(My_Year+20) + 0.5 , Obs.obs , pch=4 )  #add points of actual mortality
legend("bottomleft", legend = c("Observed", "Central forecast", "95% CI Prediction error","Actual"),
       lty = c(-1, 1, 2, 0), pch = c(1,-1,-1,4), lwd = 2, bty = "n")  #new legend

#Parameter error
Sarima.out = sarima(Kappa.hat, p=0, d=1, q=0, details = FALSE)  #sarima on kappa
Sarima.out$ttable     #test significance
SE.Param = (1:N.Ahead) * sqrt(Sarima.out$fit$sigma2/(n.y-1))  #standard error of kappa
Kappa.Up.Param = Central + Z*SE.Param   #Upper CI
Kappa.Dn.Param = Central - Z*SE.Param   #Lower CI
Forecast = Alpha.hat + Beta.hat %*% t(Central)    #fit Lee-Carter model
Forecast.Up.Param = Alpha.hat + Beta.hat %*% t(Kappa.Up.Param)  #Upper CI
Forecast.Dn.Param = Alpha.hat + Beta.hat %*% t(Kappa.Dn.Param)  #Lower CI

#Graphical output at age 65
LCForeParam <- function(){    #Lee-Carter forecast model
  Plot.Age = 65
  Plot.Row = Plot.Age - min(AGE) + 1
  par(mar=c(4.2, 4, 1, 0.5), mgp=c(3, 1, 0), las=1, cex=1)
  #plot Lee-Carter model for existing points
  plot(YEAR+0.5, Obs[Plot.Row, ], axes = FALSE, xlab = "Year", ylab = "log(mortality)",
       xlim = c(1950, 2003),ylim = c(min(Forecast.Dn.Param[Plot.Row, ]), max(Obs[Plot.Row, ]) ),
       main="Lee-Carter model fitted to males, CMI, ages 40-90, from 1950 to 1983")
  axis(1,las = 1, at = seq(1950, 2003, by = 10), tcl = -0.4)
  axis(2, seq(-6,-3, by = 0.5), tcl = -0.4)
  #add forecasted log of mortality and its confidence interval
  lines(Range+0.5, Forecast[Plot.Row, ], lwd = 2)
  lines(Range+0.5, Forecast.Up.Param[Plot.Row, ], lty = 2, lwd = 2)
  lines(Range+0.5, Forecast.Dn.Param[Plot.Row, ], lty = 2, lwd = 2)
}

LCForeParam()   #plot graph
legend("bottomleft", legend = c("Observed", "Central forecast", "95% CI Parameter error"),
       lty = c(-1, 1, 2, 2), pch = c(1,-1,-1,-1), lwd = 2, bty = "n")   #legend

LCForeParam()   #plot graph
points( (My_Year+1):(My_Year+20) + 0.5 , Obs.obs , pch=4 )    #add points of actual mortality
legend("bottomleft", legend = c("Observed", "Central forecast", "95% CI Parameter error","Actual"),
       lty = c(-1, 1, 2, 0), pch = c(1,-1,-1,4), lwd = 2, bty = "n")  #legend

Fitted.M.hat[(65-AGE[1]+1),(My_Year-YEAR[1]+1)]   #fitted value at age 65, year 1983
Forecast[(65-AGE[1]+1),ncol(Forecast)]    #projected value at age 65, year 2003



