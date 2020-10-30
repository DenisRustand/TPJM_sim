
# 1- This code shows how to simulate a dataset assuming a marginal two-part joint model
# 2- The estimation of the marginal two-part joint model is then done with frailtypack

library(frailtypack)
library(PermAlgo)
library(mvtnorm)

###########
###  1  ### Simulation of a dataset
###########
nsujet=150 # number of indivduals
# fixed effects of the model
## Binary part
alpha_0=6 # intercept
alpha_1=-3 # slope
alpha_2=1 # baseline treatment
alpha_3=-2 # treatment x slope
## Continuous part
beta_0=1.5 # intercept
beta_1=-0.5 # slope
beta_2=0.3 # baseline treatment
beta_3=0.3 # treatment x slope
sigma_e=0.3 # error term
## Survival part
gamma_1=-0.2 # treatment

gapLongi=0.25 # gap between repeated measurements of the biomarker
gap=0.001 # used to generates a lot of biomarker measurements because 
# the permutation algorithm choses among those measuremnts to define survival time

assocCL=0.08 # current-level association between two-part model and survival model

followup=4 # duration of the study

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

mestime=seq(0,followup,gap) # measurement times
timej=rep(mestime, nsujet) # time column 
nmesindiv=followup/gap+1 # number of individual measurements

nmesy= nmesindiv*nsujet# number of longi measurements
idY<-as.factor(rep(1:nsujet, each=nmesindiv)) # id

# random effects generation
MVnorm <- rmvnorm(nsujet, rep(0, 3), Sigma)
a_i = MVnorm[,3] # binary intercept
a_iY <- rep(a_i, each=nmesindiv)
b_i = MVnorm[,1] # continuous intercept
b_iY <- rep(b_i, each=nmesindiv)
bt_i = MVnorm[,2] # continuous slope
bt_iY <- rep(bt_i, each=nmesindiv)

e_ij = rnorm(nmesy,mean=0, sigma_e) # error term (continuous part)

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
# linear predictor (free from error term, for the association with survival)
linPredContTrue <- beta_0+b_iY+(beta_1+bt_iY)*timej+beta_2*trtY+beta_3*timej*trtY 
# location parameter of the lognormal distribution of positive values
mu=linPredCont-log(probaBin)-sigma_e^2/2 
muTrue=linPredCont-log(probaBin)
Y <- rlnorm(length(mu), meanlog = mu, sdlog = sigma_e) # observed positive value
Yt <- rlnorm(length(mu), meanlog = muTrue, sdlog = 0) # error-free positive value
# include zeros in the biomarker distribution
Yobs = (ifelse(B==1, Y, 0))
Ytrue = (ifelse(B==1, Yt, 0))

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
timesLongi=mestime[which(round(mestime,3) %in% round(c(seq(0,followup,by=gapLongi)),3))] # visit times
longDat <- longDat3[longDat3$timej%in%timesLongi,]
survDat$id <- as.integer(survDat$id)
longDat$id <- as.integer(longDat$id)

print(head(longDat, 20))
print(head(survDat, 20))
print(str(survDat))
print(str(longDat))
print(summary(survDat))
print(summary(longDat))

###########
###  2  ### Estimation of the marginal two-part joint model
###########
numInt=500 # number of integration points
fsurv <- Surv(deathTimes, d)~trt # survival model formula
flon <- Y~timej*trtY # continuous model formula
fbin <- Y~timej*trtY # binary model formula
# kappa value (smoothing) chosen by cross-validation
tte <- frailtyPenal(fsurv,
                    n.knots=5,kappa=0, data=survDat,cross.validation = T)
kap <- round(tte$kappa,2)

MTPJM <- longiPenal(fsurv, flon, data=survDat,data.Longi = longDat, random = c("1", "timej"),
          formula.Binary=fbin, random.Binary=c("1"),
          GLMlog=T, # logarithm link for the distribution of positive values
          MTP=T, # Trigger for marginal two-part model (set to FALSE for a conditional two-part model)
          timevar="timej",id = "id", link = "Current-level", left.censoring = F,seed.MC=1, 
          n.knots = 5, kappa = kap,hazard="Splines-per",maxit=200, 
          method.GH="Monte-carlo", n.nodes=numInt)
print(MTPJM)
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            