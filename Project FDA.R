# Load neccessary libraries
library(fda)
library(depthTools)
# Read the data from CSV file and store it in a data frame
data_frame <- read.csv("~/Documents/Uni/Functional Data/Temp.csv")
# Remove the first column (assuming it's an index column)
data_frame <- data_frame[,-1] 
df <- data_frame
# Convert the data frame to a matrix and transpose it
dataonly <- t(as.matrix(df))

time_seq <- seq(2000, 2021, 1)
# Close any open graphical devices
#dev.off()
# Plot the change in temperature data using matplot
matplot(x = time_seq, y = df, type='l', lty=1, xlab='Years', ylab='Percent of Change in C', main='% of Changes in Temprature')

# Create a B-spline basis for modeling changes in temperature trends
defbasisTray <- create.bspline.basis(c(2000, 2021), norder=11)
plot(defbasisTray)

# Evaluate the fitted temperature data at given time points for country with index of 5
time_seq <- seq(2000, 2021, 1)
t_dataonly = t(dataonly)
matrix_t_dataonly = as.matrix(t_dataonly[,5])
funct_tray <- Data2fd(argvals = time_seq, y = matrix_t_dataonly, basisobj = defbasisTray)
tray_fitted <- eval.fd(time_seq, funct_tray)

# Set up plotting parameters to display all plots in one pane
par(mfrow = c(1, 1))

# Plotting fitted data and actual data points in one pane
plot(time_seq, tray_fitted, type='l', ylim=c(-0.5, 4.5), xlab='Years', ylab='Fitted')
points(time_seq, t_dataonly[,5], col='blue', pch=20)  

# This section will present you with impact of the number of orders on fitting curve
# 
# # Define the range of norder values
# norder_values <- 3:11
# 
# # Create a for loop to iterate over each norder value
# for (norder in norder_values) {
#   # Create basis with current norder value
#   defbasisTray <- create.bspline.basis(c(2000, 2021), norder = norder)
#   
#   # Evaluate the fitted temperature data at given time points for country with index of 5
#   time_seq <- seq(2000, 2021, 1)
#   t_dataonly <- t(dataonly)
#   matrix_t_dataonly <- as.matrix(t_dataonly[,5])
#   funct_tray <- Data2fd(argvals = time_seq, y = matrix_t_dataonly, basisobj = defbasisTray)
#   tray_fitted <- eval.fd(time_seq, funct_tray)
#   
#   # Set up plotting parameters to display each plot in a separate pane
#   par(mfrow = c(1, 1))
#   
#   # Plot fitted data and actual data points in each pane
#   plot(time_seq, tray_fitted, type = 'l', ylim = c(-0.5, 4.5), xlab = 'Years', ylab = 'Fitted')
#   points(time_seq, t_dataonly[,5], col = 'blue', pch = 20)
#   
#   # Add title indicating the current value of norder
#   title(main = paste("norder =", norder))
# }


# Define a Fourier basis for modeling temperature trends
four_basis<-create.fourier.basis(rangeval = c(2000,2021), nbasis=4, period = 22)
funct_tray <- Data2fd(argvals = time_seq, y = matrix_t_dataonly, basisobj = four_basis)
tray_fitted <- eval.fd(time_seq, funct_tray)

# Plotting fitted data and actual data points
plot(time_seq, tray_fitted, type='l', ylim=c(-0.5, 4.5), xlab='Years', ylab='Fitted')
points(time_seq, matrix_t_dataonly, col='red')
title(main = "Fourier basis with nbasis = 4")


# Create a B-spline basis with lambda penalty term
norder = 6
nbasis = 22 + norder - 2
temp_basis = create.bspline.basis(c(2000,2021),
                                   nbasis, norder, time_seq)
temp_fdPar = fdPar(temp_basis, Lfdobj = 4 ,lambda = 5)
temp_smooth = smooth.basis(time_seq, matrix_t_dataonly,temp_fdPar)

#extract the estimated functional data
temp_fd=temp_smooth$fd 
# Evaluate the fitted temperature data at given time points
temp_fitted<-eval.fd(time_seq,temp_fd)
matplot(time_seq, temp_fitted, type='l', lty=1, pch='o',ylim=c(-0.5, 4.5), xlab='Years',ylab='Fitted')
points(time_seq, matrix_t_dataonly, col='purple',pch=20)
title(main = "Fourier basis with nbasis = 4, and lambda = 5")


# Calculate the second derivative (acceleration) of the smoothed temperature data
accelfdUN = deriv.fd(temp_fd, 2)
accel=eval.fd(time_seq,accelfdUN)
accelmeanfdUN = rowMeans(accel)
#join the mean to the data matrix
a=cbind(accel, accelmeanfdUN) 
# Plot the acceleration data
matplot(time_seq, a, type='l', lty=1, pch='o', xlab='Years',ylab='Acceleration/Second Derivative', col=c(rep(3,38),1), ylim=c(-2,2))
title(main = "Acceleration of the Estimated Curve")


landmark_registration <- function(lambda, norder, a,b,c) {
  # Define time sequence
  basis = norder
  time_seq <- seq(2000, 2021, len = 22)
  nbasis <- 22 + norder - 2
  temp_basis <- create.bspline.basis(c(2000, 2021), nbasis, norder, time_seq)
  temp_fdPar <- fdPar(temp_basis, Lfdobj = 4, lambda = lambda)
  temp_smooth <- smooth.basis(time_seq, t_dataonly, temp_fdPar)
  temp_fd <- temp_smooth$fd 
  temp_fitted <- eval.fd(time_seq, temp_fd)
  accelfdUN <- deriv.fd(temp_fd, 2)
  accel <- eval.fd(time_seq, accelfdUN)
  for (i in c(1:38)){
    landmax[i]=which.max(accel[c(a:b),i])+c
    landmin[i]=which.min(accel[c(a:b),i])+c
  }
  meanaccel <- rowMeans(accel)
  landmax0 <- which.max(meanaccel[a:b]) + c
  landmin0 <- which.min(meanaccel[a:b]) + c
  land_i <- cbind(time_seq[landmax], time_seq[landmin]) 
  land_0 <- c(time_seq[landmax0], time_seq[landmin0]) 
  wbasisLM <- create.bspline.basis(c(2000, 2021), basis, 4, c(2000 , time_seq[landmax0],time_seq[landmin0], 2021))
  WfdLM <- fd(matrix(0, basis, 1), wbasisLM)
  WfdParLM <- fdPar(WfdLM, 1, 0.5) 
  accelreg <- landmarkreg(accelfdUN, ximarks = land_i, x0marks = land_0, WfdPar = WfdParLM)
  return(accelreg)
}

result <- landmark_registration(lambda = 4, norder = 6, 2,17,3)

accelreg = result

accelregfd=eval.fd(temp_fine,accelreg$regfd)
accelmeanreg = rowMeans(accelregfd)

# Plot the data and the new mean
a=cbind(accelregfd, accelmeanreg) 
matplot(temp_fine, a, type='l', lty=1, pch='o', 
        xlab='Years',ylab='Acceleration', 
        main='Landmark Registration', col=c(rep(3,38),1), 
        ylim=c(-1.5,1))
# Plot the warping function
warpfd=accelreg$warpfd
warpval=eval.fd(temp_fine,warpfd)
plot(temp_fine,warpval[,1],type = 'l', col='blue', xlab='Years',
     ylab='h(Years)', main='Warping function')
lines(temp_fine,temp_fine,type='l',col='red')

##################
temp_smooth <- smooth.basis(time_seq, t_dataonly, temp_fdPar)
wbasisCR = create.bspline.basis(c(2000,2021), 15,4)
# Initialize with a set of null functions
Wfd0CR = fd(matrix(0,15,1),wbasisCR) 
# Use a penalty of order 1 but with lambda=1
WfdParCR = fdPar(Wfd0CR, 1, 1) 
target=mean.fd(temp_smooth$fd)
# Perform continuous registration
regList = register.fd(target,temp_smooth$fd, WfdParCR)
regfd=eval.fd(time_seq,regList$regfd)
meanreg = rowMeans(regfd)

# Plot the data and the new mean
a=cbind(regfd, meanreg) 
matplot(time_seq, a, type='l', lty=1, pch='o',
        xlab='Years',ylab='Acceleration', 
        main='Continuous Registration', col=c(rep(3,38),1))

# Get the warping function from continuous registration
warpfd1=regList$warpfd
warpval1=eval.fd(time_seq,warpfd1)
plot(time_seq,warpval1[,1],type = 'l', col='blue', xlab='Years',
     ylab='h(Years)', main='Warping function')
lines(time_seq,time_seq,type='l',col='red')


# Perform PCA on the registered temperature data
PCout=pca.fd(regList$regfd, nharm = 2)
print(PCout$varprop)
plot.pca.fd(PCout)

# Plot the harmonics of principal components
plot.fd(PCout$harmonics, col=c(4,2), 
        main='blue=PC1, red=PC2', xlab = 'Years', 
        ylab = 'Principal components')

# Get the variance-maximizing principal components
varmax_pcout=varmx.pca.fd(PCout)
plot.pca.fd(varmax_pcout)
plot.fd(varmax_pcout$harmonics, col=c(4,2), 
        main='blue=PC1, red=PC2', xlab = 'Years', 
        ylab = 'PCs after VARMAX')

# Read productivity data from CSV file
prod = read.csv("~/Documents/Uni/Functional Data/Prod.csv")
prod$Prod <- as.numeric(gsub("[[:space:]]", "", prod$Prod))
# Create a bar plot
barplot(prod$Prod, names.arg = prod$Country, col = "skyblue", las = 2,
        main = "Productivity by Country", ylab = "Productivity")
mat <- matrix(prod$Prod, ncol = 1, dimnames = list(prod$Country, NULL))

prod_mat = (mat)

## Compute the log10 of the total annual model_freg_proditation
## for each Region
annual_prod <- log10(as.numeric(prod_mat))
names(annual_prod) = c(prod$Country) 
smallbasis  <- create.fourier.basis(c(2000, 2021), 5)
# Build the functional data
tempfd <- smooth.basis(time_seq,
                       t_dataonly, smallbasis)
tempfd = tempfd$fd
templist = vector("list",2)
templist[[1]] = rep(1,38)
templist[[2]] = tempfd
conbasis = create.constant.basis(c(0,21))
betabasis = create.fourier.basis(c(0,21),5)
betalist = vector("list",2)
betalist[[1]] = conbasis
betalist[[2]] = betabasis

# Perform functional regression
model_freg_prod <- fRegress(annual_prod ~ tempfd)
annual_prod.fit3 <- model_freg_prod$yhatfdobj
# Plot the fit
plot(model_freg_prod$betaestlist[[2]], xlab="Years",
     ylab="Beta for temperature", main='Regression with input setup')
# Plot the data and the fit
plot(annual_prod.fit3, annual_prod, type="p", pch=20, col = "red", xlab='Predicted y', ylab='Observed y')
lines(annual_prod.fit3, annual_prod.fit3)


# Let's assess the quality of fit
annual_prodhat1 = model_freg_prod$yhatfdobj
annual_prodres1 = annual_prod - annual_prodhat1
SSE1.1 = sum((annual_prodres1)^2)
SSE0 = sum((annual_prod - mean(annual_prod))^2)
# We can now compute the squared multiple correlation 
# and the usual F-ratio 
RSQ1 = (SSE0-SSE1.1)/SSE0
Fratio1 = ((SSE0-SSE1.1)/4)/(SSE1.1/48)
RSQ1
Fratio1
#Compute the (approximate) pvalue for the test
# H0: all coefficients are identically 0
pvalue=1-pf(Fratio1,4,48)
pvalue

# Perform functional depth analysis
tempav=t_dataonly
time_range  <- time_seq
year_range <- c(2000,2021)
year_basis <- create.fourier.basis(year_range, 5)
smoothList <- with(data_frame, smooth.basis(time_seq,
                                    t_dataonly,
                                    year_basis, fdnames=list("Years", "EU Region", " % change Deg C")))
temp_fd<-smoothList$fd
plot(temp_fd)
b1<-boxplot(temp_fd, method = "MBD")
b2<-boxplot(temp_fd, method = "BD2")
b3<-boxplot(temp_fd, method = "Both")

b1$medcurve
b2$medcurve
b3$medcurve

# Load Regional Codes and Perform Wilcoxon test between different regions
regional = read.csv("~/Documents/Uni/Functional Data/Regional Code.csv")
mat_reg <- matrix(regional$Region, ncol = 1, dimnames = list(regional$Country, NULL))
region = mat_reg
DM<-b1$depth
index_arctic<-which(region=='Southern Europe')
D_arctic<-DM[index_arctic]
D_other<-DM[-index_arctic]
wilcox.test(D_arctic,D_other)

index_arctic<-which(region=='Western Europe')
D_arctic<-DM[index_arctic]
D_other<-DM[-index_arctic]
wilcox.test(D_arctic,D_other)

index_arctic<-which(region=='Eastern Europe')
D_arctic<-DM[index_arctic]
D_other<-DM[-index_arctic]
wilcox.test(D_arctic,D_other)

index_arctic<-which(region=='Northern Europe')
D_arctic<-DM[index_arctic]
D_other<-DM[-index_arctic]
wilcox.test(D_arctic,D_other)

index_arctic<-which(region=='Central Europe')
D_arctic<-DM[index_arctic]
D_other<-DM[-index_arctic]
wilcox.test(D_arctic,D_other)


# FANOVA
Years=c(2000:2021) 
matplot(Years,tempav,type='l', lty=1, pch='o', xlab='Years',ylab='Temperature', main='Temperature profile for t')


zones <- c("Constant","Western Europe", "Central Europe", "Southern Europe", "Northern Europe","Eastern Europe")

zlabels <- vector("list",6) 
zlabels[[1]] <- "Constant" 
zlabels[[2]] <- "Western Europe" 
zlabels[[3]] <- "Central Europe" 
zlabels[[4]] <- "Southern Europe" 
zlabels[[5]] <- "Northern Europe"
zlabels[[6]] <- "Eastern Europe"

# Add an empty variable to the dataset (useful later)
data_frame_new = data_frame
station <- data_frame_new$station

souindex <- which(region=='Southern Europe')
wesindex <- which(region=='Western Europe')
easindex <- which(region=='Eastern Europe')
norindex <- which(region=='Northern Europe')
cenindex <- which(region=='Central Europe')

# Set up a design matrix having a column for #the grand mean, and a column for each #climate zone effect.

# Add a dummy constraint observation
zmat <- matrix(0,38,6) 
zmat[        ,1] <- 1 
zmat[souindex,2] <- 1 
zmat[wesindex,3] <- 1 
zmat[easindex,4] <- 1 
zmat[norindex,5] <- 1
zmat[cenindex,6] <- 1


#attach a row of 0, 1, 1, 1, 1 to define the constraint
z_const <- matrix(1,1,6) 
z_const[1] <- 0 
zmat <- rbind(zmat, z_const) 
zmat

# Predicting Temperature from regional zone
smoothList <- with(data_frame_new, smooth.basis(time_seq,
                                        t_dataonly,
                                        year_basis, fdnames=list("Years", "EU Region", " % change Deg C")))
#vector with a breakpoint in the middle of each Years: created for graphical reasons
time_range  <- (2000:2021) + 0.5
year_range <- c(2000,2021)
year_basis <- create.fourier.basis(year_range, 15)
smoothList <- with(data_frame_new, smooth.basis(time_range,
                                        t_dataonly,
                                        year_basis, fdnames=list("Years", "EU Region", " % change Deg C")))
# Extract the built functional data and the matrix y2cMap which maps the data to the coefficients
Yearstempfd <- smoothList$fd 
tempy2cMap <- smoothList$y2cMap 
Yearstempfd$fdnames <- list(NULL, station, NULL)

coef <- Yearstempfd$coefs 
coef36 <- cbind(coef,matrix(0,15,1)) 
Yearstempfd$coefs <- coef36

p <- 6 
xfdlist <- vector("list",p)
for (j in 1:p) xfdlist[[j]] <- zmat[,j]

# Set up the basis for the regression functions
nbetabasis <- 15 
betabasis <- create.fourier.basis(year_range, nbetabasis)
harmaccelLfd = vec2Lfd(c((2*pi/22)^2,0,1), year_range)
betafd <- fd(matrix(0,nbetabasis,1), betabasis) 
estimate <- TRUE 
lambda <- 0 
betafdPar <- fdPar(betafd, harmaccelLfd, lambda, estimate)
betalist <- vector("list",p) 
for (j in 1:p) betalist[[j]] <- betafdPar
fRegressList <- fRegress(Yearstempfd, xfdlist, betalist)
betaestlist <- fRegressList$betaestlist

par(mfrow=c(2,3))
for (j in 1:p) { betaestParfdj <- betaestlist[[j]]
  plot(betaestParfdj$fd, xlab="Years", ylab="% change Temperature (deg C)")
  title(zlabels[[j]])
}

yhatfdobj <- fRegressList$yhatfdobj
time_range  <- (2000:2021)

yhatmat <- eval.fd(time_range, yhatfdobj) 
ymat <- eval.fd(time_range, Yearstempfd) 
tempresmat <- ymat[,1:39] - yhatmat[,1:39] 
SigmaE <- var(t(tempresmat))

# Plot estimated coefficients with a 2-standard errors confidence band
par(mfrow=c(1,1)) 
stddevE <- sqrt(diag(SigmaE)) 
plot(time_range, stddevE, type="l", xlab="Years", ylab="Standard error (% change in deg C)")
stderrList <- fRegress.stderr(fRegressList, tempy2cMap, SigmaE)
betastderrlist <- stderrList$betastderrlist

# Plot regression function standard errors
par(mfrow=c(2,3)) 
for (j in 1:p) {betastderrj <- eval.fd(time_range, betastderrlist[[j]]) 
plot(time_range, 
     type="l",lty=1, xlab="Years", ylab="Reg. Coeff.",
     betastderrj,main=zones[j])

}
par(mfrow=c(1,1)) 
for (j in 1:p) {
  betafdParj <- betaestlist[[j]] 
  betafdj <- betafdParj$fd 
  betaj <- eval.fd(time_range, betafdj) 
  betastderrj <- eval.fd(time_range, betastderrlist[[j]]) 
  matplot(time_range, cbind(betaj, betaj+2*betastderrj, betaj-2*betastderrj),
          type="l",lty=c(1,4,4), xlab="Years", ylab="Reg. Coeff.",
          main=zones[j]
  )
}
# Plot temperature residuals
lambda <- 22
fdParobj <- fdPar(year_basis, harmaccelLfd, lambda)
smoothList <- smooth.basis(time_range, tempresmat, fdParobj)
tempresfdobj <- smoothList$fd
par(mfrow=c(1,1)) 
plot(tempresfdobj)
