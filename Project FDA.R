df <- read.csv("/Users/amedodavoodi/Desktop/Temp.csv")
df <- df[,-1]  # Check if removing the first column is intentional
au12 <- df
dataonly <- t(as.matrix(au12))

matplot(dataonly, type='l', lty=1, xlab='time', ylab='Percent of Change in Temprature', main='Trends of changes in Temprature - EU')

library(fda)

# Define the basis for tray47
defbasisTray <- create.bspline.basis(c(2000, 2021), norder=4)
plot(defbasisTray)

years <- seq(2000, 2021, 1)

t_dataonly = t(dataonly)
x = as.matrix(t_dataonly[,10])

funct_tray47 <- Data2fd(argvals = years, y = x, basisobj = defbasisTray)
tray47_fitted1 <- eval.fd(years, funct_tray47)

# Plotting fitted data and actual data points
plot(years, tray47_fitted1, type='l', ylim=c(-0.5, 4.5), xlab='t', ylab='Tray 47 level')
points(years, x, col='lightblue')

def_fourier_b<-create.fourier.basis(rangeval = c(2000,2021), nbasis=4, period = 22) 


funct_tray47 <- Data2fd(argvals = years, y = x, basisobj = def_fourier_b)
tray47_fitted1 <- eval.fd(years, funct_tray47)

# Plotting fitted data and actual data points
plot(years, tray47_fitted1[1:22], type='l', ylim=c(-0.5, 4.5), xlab='t', ylab='Tray 47 level')
points(years, t_dataonly[1:22], col='lightblue')



##########################
#years
#x > for one var
#dataonly = all vars

norder = 6
nbasis = 22 + norder - 2
heightbasis = create.bspline.basis(c(2000,2021),
                                   nbasis, norder, years)
heightfdPar = fdPar(heightbasis, Lfdobj = 4 ,lambda = 0.5)
height_smooth = smooth.basis(years, t_dataonly,heightfdPar)
heightfd=height_smooth$fd #extract the estimated functional data
height_fitted<-eval.fd(years,heightfd)
matplot(years, height_fitted, type='l', lty=1, pch='o', xlab='Years',ylab='Temp Change')

##################


accelfdUN = deriv.fd(heightfd, 2)
accel=eval.fd(years,accelfdUN)
accelmeanfdUN = rowMeans(accel)
a=cbind(accel, accelmeanfdUN) #join the mean to the data matrix
# for the plot
matplot(years, a, type='l', lty=1, pch='o', xlab='years',ylab='acceleration', col=c(rep(3,38),1), ylim=c(-2,2))


#################
# 2005 - 2015

agefine= seq(2000,2021,len=22)
accel=eval.fd(agefine,accelfdUN)
#agefine[6]=2005 and agefine[16]=2015: components
#will be our time extremes of interest to identify landmarks

#initialize vectors for the location of maxima and minima
landmax=rep(0,38) 
landmin=rep(0,38)
for (i in c(1:38)){
  landmax[i]=which.max(accel[c(6:16),i])+5
  landmin[i]=which.min(accel[c(6:16),i])+5
}

meanaccel=rowMeans(accel)
landmax0=which.max(meanaccel[c(6:16)])
landmin0=which.min(meanaccel[c(6:16)])


land_i=cbind(agefine[landmax],agefine[landmin]) #matrix of landmarks of data
land_0=c(agefine[landmax0],agefine[landmin0]) #landmarks of target



#create a basis for warping
basss = 6
wbasisLM=create.bspline.basis(c(2000,2021),basss, 4, c(2000,agefine[landmin0],agefine[landmax0],2021))

#to be used as initial guess in the landmark registration
WfdLM = fd(matrix(0,basss,1),wbasisLM)
#smoothing with a penalty of order 1
WfdParLM = fdPar(WfdLM,1,0.5) 

#apply landmarkreg
accelreg=landmarkreg(accelfdUN,ximarks=land_i,x0marks=land_0,WfdPar=WfdParLM)

accelregfd=eval.fd(agefine,accelreg$regfd)
accelmeanreg = rowMeans(accelregfd)
#plot the data and the new mean
a=cbind(accelregfd, accelmeanreg) 
matplot(agefine, a, type='l', lty=1, pch='o', 
        xlab='age',ylab='acceleration', 
        main='landmark registration', col=c(rep(3,38),1), 
        ylim=c(-1,1))

warpfd=accelreg$warpfd
warpval=eval.fd(agefine,warpfd)
plot(agefine,warpval[,1],type = 'l', col='blue', xlab='age',
     ylab='h(age)', main='warping function')
lines(agefine,agefine,type='l',col='grey')



##################

wbasisCR = create.bspline.basis(c(2000,2021),15, 5)
Wfd0CR = fd(matrix(0,15,1),wbasisCR) #initialize with a set of null functions
WfdParCR = fdPar(Wfd0CR, 1, 1) #use a penalty of order 1 but with lambda=1
target=mean.fd(height_smooth$fd)
regList = register.fd(target,height_smooth$fd, WfdParCR)
regfd=eval.fd(years,regList$regfd)
meanreg = rowMeans(regfd)
#plot the data and the new mean
a=cbind(regfd, meanreg) 
matplot(years, a, type='l', lty=1, pch='o',
        xlab='age',ylab='registered height', 
        main='continuous registration', col=c(rep(3,38),1))
##################

warpfd1=regList$warpfd
warpval1=eval.fd(years,warpfd1)
plot(years,warpval1[,1],type = 'l', col='blue', xlab='age',
     ylab='h(age)', main='warping function')
lines(years,years,type='l',col='grey')
#######################

PCout=pca.fd(regList$regfd, nharm = 2)
print(PCout$varprop)
plot.pca.fd(PCout)


#######################

plot.fd(PCout$harmonics, col=c(4,2), 
        main='blue=PC1, red=PC2', xlab = 'age', 
        ylab = 'Principal components')
########################

varmax_pcout=varmx.pca.fd(PCout)
plot.pca.fd(varmax_pcout)

plot.fd(varmax_pcout$harmonics, col=c(4,2), 
        main='blue=PC1, red=PC2', xlab = 'age', 
        ylab = 'PCs after VARMAX')
############################

prod = read.csv("/Users/amedodavoodi/Desktop/Prod.csv")
prod$Prod <- as.numeric(gsub("[[:space:]]", "", prod$Prod))


# Create a bar plot
barplot(prod$Prod, names.arg = prod$Country, col = "skyblue", las = 2,
        main = "Productivity by Country", ylab = "Productivity")


mat <- matrix(prod$Prod, ncol = 1, dimnames = list(prod$Country, NULL))

prod_mat = (mat)

##
## Compute the log10 of the total annual precipitation
## for each station
##
annualprec <- log10(as.numeric(prod_mat))
names(annualprec) = c(prod$Country) 


smallbasis  <- create.fourier.basis(c(2000, 2021), 5)
# Build the functional data
tempfd <- smooth.basis(years,
                       t_dataonly, smallbasis)
tempfd = tempfd$fd

templist = vector("list",2)
templist[[1]] = rep(1,38)
templist[[2]] = tempfd

##

conbasis = create.constant.basis(c(0,21))
betabasis = create.fourier.basis(c(0,21),5)

betalist = vector("list",2)
betalist[[1]] = conbasis
betalist[[2]] = betabasis


precip <- fRegress(annualprec ~ tempfd)

annualprec.fit3 <- precip$yhatfdobj


#plot the fit
plot(precip$betaestlist[[2]], xlab="Day",
     ylab="Beta for temperature", main='Regression with input setup')
#  plot the data and the fit
plot(annualprec.fit3, annualprec, type="p", pch="o", xlab='predicted y', ylab='observed y')
lines(annualprec.fit3, annualprec.fit3)

##################################################

#Let's assess the quality of fit

annualprechat1 = precip$yhatfdobj
annualprecres1 = annualprec - annualprechat1
SSE1.1 = sum((annualprecres1)^2)
SSE0 = sum((annualprec - mean(annualprec))^2)

#We can now compute the squared multiple correlation 
#and the usual F-ratio 
RSQ1 = (SSE0-SSE1.1)/SSE0
Fratio1 = ((SSE0-SSE1.1)/4)/(SSE1.1/48)
RSQ1
Fratio1

#Compute the (approximate) pvalue for the test
#H0: all coefficients are identically 0
pvalue=1-pf(Fratio1,df1=4,df2=48)
pvalue


###############################

library(fda)
library(depthTools)

tempav=t_dataonly
daytime  <- years
dayrange <- c(2000,2021)
daybasis15 <- create.fourier.basis(dayrange, 5)
smoothList <- with(df, smooth.basis(years,
                                                 t_dataonly,
                                                 daybasis15, fdnames=list("Day", "Station", "Deg C")))
temp_fd<-smoothList$fd
plot(temp_fd)


b1<-boxplot(temp_fd, method = "MBD")

b2<-boxplot(temp_fd, method = "BD2")
b3<-boxplot(temp_fd, method = "Both")

b1$medcurve
b2$medcurve
b3$medcurve


regional = read.csv("/Users/amedodavoodi/Desktop/Regional Code.csv")
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


##############################
# FANOVA

day=c(2000:2021) 
matplot(day,tempav,type='l', lty=1, pch='o', xlab='day',ylab='temperature', main='Temperature profile for t')


zones <- c("Constant","Western Europe", "Central Europe", "Southern Europe", "Northern Europe","Eastern Europe")

zlabels <- vector("list",6) 
zlabels[[1]] <- "Constant" 
zlabels[[2]] <- "Western Europe" 
zlabels[[3]] <- "Central Europe" 
zlabels[[4]] <- "Southern Europe" 
zlabels[[5]] <- "Northern Europe"
zlabels[[6]] <- "Eastern Europe"

#add an empty variable to the dataset (useful later)
df_new = df
station <- df_new$station

souindex <- which(region=='Southern Europe')
wesindex <- which(region=='Western Europe')
easindex <- which(region=='Eastern Europe')
norindex <- which(region=='Northern Europe')
cenindex <- which(region=='Central Europe')

# Set up a design matrix having a column for #the grand mean, and a column for each #climate zone effect.

#Add a dummy constraint observation

zmat <- matrix(0,38,6) 
zmat[ ,1] <- 1 
zmat[souindex,2] <- 1 
zmat[wesindex,3] <- 1 
zmat[easindex,4] <- 1 
zmat[norindex,5] <- 1
zmat[cenindex,6] <- 1


#attach a row of 0, 1, 1, 1, 1 to define the constraint

z36 <- matrix(1,1,6) 
z36[1] <- 0 
zmat <- rbind(zmat, z36) 
zmat

#######
smoothList <- with(df, smooth.basis(years,
                                    t_dataonly,
                                    daybasis15, fdnames=list("Day", "Station", "Deg C")))

daytempfd <- smoothList$fd 
tempy2cMap <- smoothList$y2cMap 
daytempfd$fdnames <- list(NULL, station, NULL)

coef <- daytempfd$coefs 
coef36 <- cbind(coef,matrix(0,15,1)) 
daytempfd$coefs <- coef36

p <- 6 
xfdlist <- vector("list",p)
for (j in 1:p) xfdlist[[j]] <- zmat[,j]


nbetabasis <- 13 
betabasis <- create.fourier.basis(dayrange, nbetabasis)

harmaccelLfd = vec2Lfd(c((2*pi/365)^2,0,1), dayrange)

betafd <- fd(matrix(0,nbetabasis,1), betabasis) 
estimate <- TRUE 
lambda <- 0 
betafdPar <- fdPar(betafd, harmaccelLfd, lambda, estimate)

betalist <- vector("list",p) 

for (j in 1:p) betalist[[j]] <- betafdPar

fRegressList <- fRegress(daytempfd, xfdlist, betalist)
betaestlist <- fRegressList$betaestlist

par(mfrow=c(2,3))
for (j in 1:p) { betaestParfdj <- betaestlist[[j]]

plot(betaestParfdj$fd, xlab="Day", ylab="Temperature (deg C)")

title(zlabels[[j]])

}


yhatfdobj <- fRegressList$yhatfdobj

daytime <- (0:21) + 0.5

yhatmat <- eval.fd(daytime, yhatfdobj) 
ymat <- eval.fd(daytime, daytempfd) 
tempresmat <- ymat[,1:39] - yhatmat[,1:39] 
SigmaE <- var(t(tempresmat))

par(mfrow=c(1,1)) 
stddevE <- sqrt(diag(SigmaE)) 
plot(daytime, stddevE, type="l", xlab="Day", ylab="Standard error (deg C)")


stderrList <- fRegress.stderr(fRegressList, tempy2cMap, SigmaE)

betastderrlist <- stderrList$betastderrlist

#plot regression function standard errors

par(mfrow=c(2,3)) 
for (j in 1:p) {betastderrj <- eval.fd(daytime, betastderrlist[[j]]) 
                plot(daytime, 
                     type="l",lty=1, xlab="Day", ylab="Reg. Coeff.",
                     betastderrj,main=zones[j])
  
}



par(mfrow=c(1,1)) 
for (j in 1:p) {
  betafdParj <- betaestlist[[j]] 
  betafdj <- betafdParj$fd 
  betaj <- eval.fd(daytime, betafdj) 
  betastderrj <- eval.fd(daytime, betastderrlist[[j]]) 
  matplot(daytime, cbind(betaj, betaj+2*betastderrj, betaj-2*betastderrj),
          type="l",lty=c(1,4,4), xlab="Day", ylab="Reg. Coeff.",
          main=zones[j]
          )
}


lambda <- 22
fdParobj <- fdPar(daybasis15, harmaccelLfd, lambda)

smoothList <- smooth.basis(daytime, tempresmat, fdParobj)

tempresfdobj <- smoothList$fd

#plot temperature residuals

par(mfrow=c(1,1)) 
plot(tempresfdobj)



