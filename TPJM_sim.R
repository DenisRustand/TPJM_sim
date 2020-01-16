
# The two-part model is implemented in the last version 
# of the frailtypack package on the CRAN repository (3.1.0).

# The following code is decomposed in 3 parts:
# 1. Simulation of a dataset assuming a conditional two-part joint model 
#    (with current-level association structure)
# 2. Estimation of the true model as long as a joint naive and left-censored model
# 3. Plot of the conditional survival taking into account the effect of treatment   
#    on the risk of event through the biomarker in addition to the direct treatment 
#    effect on the risk of event.

suppressMessages(library(frailtypack))
suppressMessages(library(PermAlgo)) # for generation of death times and censoring times
suppressMessages(library(mvtnorm)) # multivariate normal generation

# Method for the numerical approximation of the integral over random-effects
methodInt="Monte-carlo" # ("Standard" for standard gauss-hermite quadrature)
# association structure for estimation ("Current-level" / "Two-part" / "Random-effects")
assoc="Current-level" 
MAXITER=100 # max iterations for Marquardt algorithm
estim_TPJM <- T # Two-part joint model estimation
estim_JMn <- F # Naive joint model estimation
estim_JMlc <- F # Left-censoring joint model estimation
plotsRes=F # plot results (survival conditional on treatment)
set.seed(1) # seed for data generation
seed_MC=1 # seed for Monte-carlo integration method

#####
# 1 # data simulation
#####
# Need to sync random number generator because of some changes in R 3.6
if("Rounding"%in%RNGkind() | "Rejection"%in%RNGkind()){
  suppressWarnings(RNGkind(sample.kind = "Rounding"))
}

nsujet=150 # number of indivduals
numInt=500 # number of integration points

# binary part fixed effects
alpha_0=6 # intercept
alpha_1=-4 # slope
alpha_2=-1 # baseline treatment
alpha_3=1 # treatment X time

# continuous part fixed effects
beta_0=4 # intercept
beta_1=-0.5 # slope
beta_2=1 # baseline treatment
beta_3=1 # treatment X time

# survival part
gamma_1=0.3 # fixed effect of treatment on the risk of event

sigma_e=0.5 # error term (gaussian)


gapLongi=0.07 # gap between longi measurements
gap=0.001 # used to generate a lot of time points because the permutation
# algorithm choses among those time points to define survival times

assocCL=0.30 # current-level association between two-part model and survival model

kno=5 #knots for splines / baseline hazard

followup=4 # follow-up time

# random effects variance and covariance
sigma_b=sqrt(0.5625) # continuous intercept
sigma_bt=sqrt(0.5625) # continuous slope
sigma_a=sqrt(4) # binary intercept
cor_bbt=-0.2 # correlation continuous intercept X slope
cor_ba=0.2 # correlation continuous intercept X binary intercept
cor_bta=0.7 # correlation continuous slope X binary intercept
cov_bbt <- sigma_b*sigma_bt*cor_bbt
cov_ba <- sigma_b*sigma_a*cor_ba
cov_bta <- sigma_bt*sigma_a*cor_bta

Sigma=matrix(c(sigma_b^2,cov_bbt,cov_ba,
               cov_bbt,sigma_bt^2,cov_bta,
               cov_ba,cov_bta,sigma_a^2),ncol=3,nrow=3)

fsurv <- Surv(deathTimes, d)~trt # survival model formula
flon <- Y~timej*trtY # continuous model formula
fbin <- Y~timej*trtY # binary model formula

mestime=seq(0,followup,gap) # measurement times
timej=rep(mestime, nsujet) # time column 
nmesindiv=followup/gap+1 # number of individual measurements

nmesy= nmesindiv*nsujet# number of longi measurements
idY<-as.factor(rep(1:nsujet, each=nmesindiv)) # id

### begin data generation
# random effects generation
MVnorm <- mvtnorm::rmvnorm(nsujet, rep(0, 3), Sigma)

a_i = MVnorm[,3] # binary intercept
a_iY <- rep(a_i, each=nmesindiv)

b_i = MVnorm[,1] # continuous intercept
b_iY <- rep(b_i, each=nmesindiv)
bt_i = MVnorm[,2] # continuous slope
bt_iY <- rep(bt_i, each=nmesindiv)

e_ij = rnorm(nmesy,mean=0, sigma_e) # error

trt=rbinom(nsujet,1, 0.5) # treatment covariate
trtY=rep(trt, each=nmesindiv)


## binary part generation
# linear predictor (binary)
linPredBin <- alpha_0+a_iY+alpha_1*timej+alpha_2*trtY+alpha_3*timej*trtY 
probaBin <- exp(linPredBin)/(1+exp(linPredBin)) # proba of zero
B <- rbinom(nmesy,1, probaBin) # zero values (binomial)


## generation of longitudinal measurements of outcome
# linear predictor (continuous)
linPredCont <- beta_0+b_iY+(beta_1+bt_iY)*timej+beta_2*trtY+beta_3*timej*trtY+e_ij 
# linear predictor (free from error term, for the association)
linPredContTrue <- beta_0+b_iY+(beta_1+bt_iY)*timej+beta_2*trtY+beta_3*timej*trtY 

# in case of negative generated continuous measurements (rarely happening)
linPredContP <- ifelse(linPredCont<(0), 0, linPredCont) 
linPredContTrueP <- ifelse(linPredCont<(0), 0, linPredContTrue) 

# include zeros in the biomarker distribution
Yobs = (ifelse(B==1, linPredContP, 0))
Ytrue = (ifelse(B==1, linPredContTrueP, 0))

#longitudinal dataset
id <- as.integer(idY)
longDat <- data.frame(id, timej, trtY, Yobs, B, Ytrue)

# longi measurements to generate survival times with permutation algorithm
matX=matrix(ncol=3, nrow=nsujet*nmesindiv)  
# treatment covariate (to evaluate treatment effect on the risk of event)
matX[,1] <- longDat[,"trtY"]
# true value of the biomarker (to evaluate effect of the biomarker on the risk of event)
matX[,2] <- longDat[,"Ytrue"] 
# observed value of the biomarker
matX[,3] <- longDat[,"Yobs"] 
eventRandom <- round(rexp(nsujet, 0.0012)+1,0) # ~80% death
censorRandom=runif(nsujet,1,nmesindiv) # uniform random censoring
Ttemp <- permalgorithm(nsujet,nmesindiv,Xmat=matX,eventRandom = eventRandom,
                       censorRandom=censorRandom,XmatNames=c("trtY", "Ytrue", "Yobs"), 
                       betas=c(gamma_1,assocCL, 0) )

# extract last line of each individual (= death/censoring time)
ligne=NULL
for(i in 1:(dim(Ttemp)[1]-1)){
  if(Ttemp[i,"Id"]!=Ttemp[i+1,"Id"]) ligne <- c(ligne, i)  
}
ligne <-c(ligne, dim(Ttemp)[1])


Ttemp2=Ttemp[ligne, c("Id","Event","Stop", "trtY")] # one line per individual
Ttemp2$deathTimes <- mestime[Ttemp2$Stop+1] # deathtimes
survDat <- Ttemp2[, c("Id", "deathTimes", "Event", "trtY")] # survival dataset
names(survDat) <- c("id", "deathTimes", "d", "trt")

longDat2 <- Ttemp[,c("Id", "Start", "trtY", "Yobs")]
longDat2$timej <- mestime[longDat2$Start+1] # measurements times of the biomarker
longDat3 <- longDat2[, c("Id", "timej", "trtY", "Yobs")]
names(longDat3) <- c("id", "timej", "trtY", "Y")
timesLongi=mestime[which(mestime%%gapLongi==0)] # visit times
longDat <- longDat3[longDat3$timej%in%timesLongi,]

survDat$id <- as.integer(survDat$id)
longDat$id <- as.integer(longDat$id)

# Datasets generated are also stored in the frailtypack
#load(survDat)
#load(longDat)

### end data generation

print(head(longDat, 20))
print(head(survDat, 20))
print(str(survDat))
print(str(longDat))
print(summary(survDat))
print(summary(longDat))

#####
# 2 # Model estimation
#####

# kappa value (smoothing) chosen by cross-validation with an univariate Cox model
tte <- frailtyPenal(fsurv,n.knots=kno,kappa=0, data=survDat,cross.validation = T)
kap <- round(tte$kappa,2);kap # smoothing parameter
if(estim_TPJM){ # computation takes ~12min with an Intel i7-4790 (8 cores, 3.60 GHz)
  TPJM <- longiPenal(fsurv, flon, data=survDat, data.Longi = longDat, 
                     random = c("1","timej"), formula.Binary=fbin, 
                     random.Binary=c("1"), timevar="timej", id = "id", 
                     link = assoc, n.knots = kno, kappa = kap,
                     hazard="Splines-per", method.GH=methodInt, 
                     n.nodes=numInt, seed.MC=seed_MC);TPJM
}


## joint naive model
if(estim_JMn){
  JMn <- longiPenal(fsurv, flon, data=survDat, data.Longi = longDat, 
                    random = c("1","timej"), timevar="timej",id = "id", 
                    link = assoc, n.knots = kno, 
                    kappa = kap,hazard="Splines-per", 
                    method.GH=methodInt, n.nodes=numInt, seed.MC=seed_MC);JMn
}


## joint left-censored model
if(estim_JMlc){
  # censoring threshold (just below smallest observed positive value)
  TRE <- min(longDat[longDat$Y!=min(longDat$Y),"Y"])-
    (min(longDat[longDat$Y!=min(longDat$Y),"Y"])/1000)
  JMlc <- longiPenal(fsurv, flon, data=survDat, data.Longi = longDat, 
                     random = c("1","timej"), timevar="timej", id = "id", 
                     link = assoc, left.censoring = TRE, n.knots = kno, 
                     kappa = kap,hazard="Splines-per",
                     method.GH=methodInt, n.nodes=numInt, seed.MC=seed_MC);JMlc
}



#####
# 3 # plots results (only Two-part model /  conditional on treatment arm)
#####

if(plotsRes){
  # Plot conditional survival from a model estimated with frailtypack
  # M-splines for the baseline hazard risk and
  # I-splines for the baseline survival (Ispline=integral(Msplines))
  # We estimate baseline survival with a numerical approximation
  
  
  # load models as R objects
  #load("~/TPJM.RData")
  TP=TPJM
  
  
  #--------------mspline-----------------------#
  #' this function generates M_i
  #' @param x time x for estimation
  #' @param tp timepoint of length n+k
  #' @param n.knot total number of knots
  #' @param k order of the spline function
  mspline = function(x,tp,n.knot,k=4){
    if(k==1){
      region = cbind(tp[1:(length(tp)-1)],tp[2:length(tp)])
      bool = as.integer(x>region[,1] & x<region[,2])
      return((1/diff(tp))[as.logical(bool)]*bool)
    }
    else{
      n=length(tp)-k
      region = cbind(tp[1:n],tp[(k+1):(k+n)])
      bool = I(x>region[,1] & x<region[,2])
      M=k*((x-tp[1:n])*(mspline(x,tp,n.knot,k-1)[1:n])+
             (tp[(k+1):(k+n)]-x)*(mspline(x,tp,n.knot,k-1)[2:(n+1)]))/
        ((k-1)*(tp[(k+1):(k+n)]-tp[1:n]))
      M_final = rep(0,n)
      M_final[bool]=M[bool]
      return(M_final)
    }
  }
  
  tpoints=seq(0,max(TP$xD),len=1000) # time points for splines estimation and biomarker values
  BH=TP$b[1:7]^2 # parameters associated to splines (n.knots+2)
  M_i=apply(as.matrix(tpoints), 1,mspline,tp=TP$zi,n.knot = 5,k=4) # M-splines
  hazardEst=t(M_i)%*%as.matrix(BH) # baseline hazard risk
  
  weights=rep(tpoints[2]-tpoints[1],len=length(tpoints)) # for the integral approximation
  hCUM=cumsum(hazardEst)*weights # baseline cumulative risk
  survEst <- exp(-hCUM) # baseline survival
  
  # biomarker value
  coefTP <- TP$coef # model parameters
  res <- NULL
  TwoPart <- function(t,trt){
    BinLinPred <- coefTP[6]+coefTP[7]*t+coefTP[8]*trt+coefTP[9]*t*trt
    ConLinPred <- coefTP[2]+coefTP[3]*t+coefTP[4]*trt+coefTP[5]*t*trt
    Prob <- exp(BinLinPred)/(1+exp(BinLinPred))
    res <- Prob*ConLinPred
    return(res)
  }
  
  survx <- tpoints # abscissas
  survy <- survEst^exp(TwoPart(survx,0)*TP$eta) # survival (arm A)
  survytrt <- survEst^exp(coefTP[1]+TwoPart(survx,1)*TP$eta) # (arm B)
  
  # confidence intervals (Monte-carlo method)
  nloop=1000 # number of Monte-Carlo curves
  Hess <- TP$varHIHtotal # Hessian matrix (splines for the baseline hazard included)
  # generation of the random points
  isCoef <- mvtnorm::rmvnorm(nloop, TP$b, Hess) 
  mc_BH=isCoef[,1:7]^2 # parameters for the splines (survival)
  M_i=apply(as.matrix(tpoints), 1,mspline,tp=TP$zi,n.knot = 5,k=4) # M-splines
  mc_hazardEst=apply(mc_BH,1,function(x) t(M_i)%*%as.matrix(x))# baseline hazard
  mc_hCUM=apply(mc_hazardEst,2,cumsum)
  
  mc_hCUMfinal = apply(mc_hCUM,2, function(x) x*weights)
  mc_survEst <- exp(-mc_hCUMfinal) # baseline survival for all the Monte-carlo curves
  
  # biomarker
  survMC=NULL
  survMCtrt=NULL
  for(i in 1:nloop){
    curve_i <- mc_survEst[,i]^exp(TwoPart(survx,0)*isCoef[i,8])
    curve_itrt <- mc_survEst[,i]^exp(isCoef[i,16]+TwoPart(survx,1)*isCoef[i,8])
    survMC <- cbind(survMC, curve_i) # survival (arm A)
    survMCtrt <- cbind(survMCtrt, curve_itrt) # survival (arm B)
  }
  
  # quantiles
  QL <- function(x) quantile(x,prob=0.025)
  QU <- function(x) quantile(x,prob=0.975)
  SCL <- apply(survMC,1,QL) # ref lower
  SCU <- apply(survMC,1,QU) # ref upper
  SCLtrt <- apply(survMCtrt,1,QL) # trt lower
  SCUtrt <- apply(survMCtrt,1,QU) # trt upper
  
  # plot
  par(mfrow=c(1,1))
  plot(survx,survy,lwd=2,xlab="time",ylab="survival",ylim=c(0,1),type='l',
       main="Survival conditional on treatment arm (TPJM current-level)")
  lines(survx,survytrt,col='red',lwd=2,lty=2)
  lines(survx,SCL)
  lines(survx,SCU)
  lines(survx,SCLtrt,col='red',lty=2)
  lines(survx,SCUtrt,col='red',lty=2)
  legend("topright",title = "Treatment", c("arm A","arm B"), lty=c(1,1),
         lwd=c(2,2),col=c("black","red"))
}




