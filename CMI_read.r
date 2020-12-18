#
#   File: CMI_read.r
#
#   CMI data
#
#   Deaths
#
x <- scan(file="CMI_Deaths.csv", what = "character", skip = 4, sep = ',')
x.mat <- matrix(x, ncol = 60, byrow = TRUE)
Dth <- x.mat[2:91, 2:60]
Dth <- matrix(as.numeric(Dth), nrow = 90, ncol = 59)
Dth.NA <- is.na(Dth)
sum(Dth[Dth.NA == TRUE])
#
Age <- 11:100
Year <- 1947:2005
dimnames(Dth) = list(Age, Year)
#
#   Exposures
#
x <- scan(file="CMI_Exposures.csv", what = "character", skip = 4, sep = ',')
x.mat <- matrix(x, ncol = 60, byrow = TRUE)
Exp <- x.mat[2:91, 2:60]
Exp <- matrix(as.numeric(Exp), nrow = 90, ncol = 59)
#
dimnames(Exp) = list(Age, Year)
#
#   Convert to central exposure
#
Exp <- Exp - Dth/2
#
#   Select ages 20-90, years 1950:2005
#
Age <- Age[10:80]
Year <- Year[4:59]
Dth <- Dth[10:80, 4:59]
Exp <- Exp[10:80, 4:59]
#
#    Tidy up
#
rm(x, x.mat, Dth.NA)

