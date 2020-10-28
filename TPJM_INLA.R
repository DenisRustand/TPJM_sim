
# 1- This code shows how to simulate a dataset assuming a conditional two-part joint model
# 2- The estimation of the conditional two-part joint model is then done with INLA

set.seed(1)
library(INLA)
inla.setOption(mkl=TRUE)

###########
###  1  ### Simulation of a dataset
###########

library(mvtnorm) # for multivariate normal generation (random-effects)
nsujet=200 # number of individuals
#binary part  
alpha_0=4 # Intercept
alpha_1=-0.5 # slope
alpha_2=-0.5 # treatment
alpha_3=0.5 # treatment x time
#continuous part
beta_0=2 # Intercept
beta_1=-0.3 # slope
beta_2=-0.3 # treatment
beta_3=0.3 # treatment x time
sigma_e=0.3 # error term (standard error)
gamma_1=0.2 # treatmentt effect on survival
# Shared random effects association between the two-part model for the biomarker and survival
phi_a=1 # random intercept (binary)
phi_b=1 # random intercept (continuous)
phi_bt=1 # random slope (continuous)
# baseline hazard scale (to generate exponential death times)
baseScale=0.2
gap=0.4# gap between longitudinal repeated measurements
followup=4 # study duration
# correlated random-effects
sigma_a=1 # random intercept (binary)
sigma_b=0.5 # random intercept (continuous)
sigma_bt=0.5 # random slope (continuous)
cor_ba=0.5 # correlation intercept (binary)/intercept (continuous)
cor_bta=0.5 # correlation intercept (binary)/slope (continuous)
cor_bbt=-0.2 # correlation continuous intercept/slope
cov_ba <- sigma_b*sigma_a*cor_ba # covariance
cov_bta <- sigma_bt*sigma_a*cor_bta
cov_bbt <- sigma_b*sigma_bt*cor_bbt
Sigma=matrix(c(sigma_a^2,cov_ba,cov_bta, # variance-covariance matrix
               cov_ba,sigma_b^2,cov_bbt,
               cov_bta,cov_bbt,sigma_bt^2),ncol=3,nrow=3)
mestime=seq(0,followup,gap) # measurement times
timej=rep(mestime, nsujet) # time column 
nmesindiv=followup/gap+1 # number of individual measurements
nmesy= nmesindiv*nsujet # number of longi measurements
id<-as.factor(rep(1:nsujet, each=nmesindiv)) # patient id
# random effects generation
MVnorm <- mvtnorm::rmvnorm(nsujet, rep(0, 3), Sigma)
a_i = MVnorm[,1] # binary intercept
a_iY <- rep(a_i, each=nmesindiv) # binary intercept (repeated for longi dataset)
b_i = MVnorm[,2] # continuous intercept
b_iY <- rep(b_i, each=nmesindiv)
bt_i = MVnorm[,3] # continuous slope
bt_iY <- rep(bt_i, each=nmesindiv)

treated <- sample(1:nsujet, nsujet/2, replace=F)
treatedFull <- NULL
for(i in 1:nsujet){
  treatedFull <- c(treatedFull, ifelse(i%in%treated, 1, 0))
}
trt= treatedFull# treatment covariate
trtY=rep(trt, each=nmesindiv)

## linear predictor (binary part)
linPredBin <- alpha_0+a_iY+alpha_1*timej+alpha_2*trtY+alpha_3*timej*trtY
probaBin <- exp(linPredBin)/(1+exp(linPredBin)) # proba of positive value
B <- rbinom(nmesy,1, probaBin) # observed zero values

## linear predictor (continuous part)
linPredCont <- beta_0+b_iY+(beta_1+bt_iY)*timej+beta_2*trtY+beta_3*timej*trtY
mu=linPredCont-sigma_e^2/2 # lognormal mean
Ypos <- rlnorm(length(mu), meanlog = mu, sdlog = sigma_e) # observed biomarker values
Y = (ifelse(B==1, Ypos, 0)) # include zeros in the biomarker distribution

## longitudinal biomarker dataset
longDat <- data.frame(id, timej, trtY, Y)

## generation of exponential death times
u <- runif(nsujet) # uniform distribution for survival times generation
deathTimes <- -(log(u) / (baseScale * exp(trt * gamma_1 + a_i*phi_a + b_i*phi_b + bt_i*phi_bt)))
d <- as.numeric(deathTimes<followup) # deathtimes indicator
## censoring individuals at end of follow-up (not at random)
deathTimes[deathTimes>=followup]=followup 
ids <- as.factor(1:nsujet)
survDat <- data.frame(id=ids,deathTimes, d, trt) # survival times dataset

## removing longi measurements after death
ind <- rep(NA, nsujet*length(mestime)) 
for (i in 1:nsujet){ 
  for(j in 1:length(mestime)){
    if(longDat[(i-1)*length(mestime)+j, "timej"]<=survDat[i,"deathTimes"]) ind[(i-1)*length(mestime)+j]=1
  }
}
longDat <- longDat[!is.na(ind),]
survDat$trt <- as.factor(survDat$trt)
longDat$trtY <- as.factor(longDat$trtY)
longDat$id <- as.integer(longDat$id)
survDat$id <- as.integer(survDat$id)
## Summary of the longitudinal and survival datasets
print(summary(survDat))
print(summary(longDat))

###########
###  1  ### Estimation of a conditional two-part joint model with R-INLA
###########

# create dataset with positive values only for the continuous part
longDatlog <- longDat[longDat$Y>0,]
nB <- length(longDat$Y) # length of binary part
nC <- length(longDatlog$Y) # length of continuous part
ns=dim(survDat)[1] # number of individuals

longDat$B <- ifelse(longDat$Y==0,0,1) # zero value indicator (binary part outcome)
longDatlog$sld <- longDatlog$Y # positive values only (continuous part outcome)
yy <- matrix(NA, ncol = 2, nrow = nB+nC)
yy[1:nB,1] <- longDat$B # binary outcome
yy[nB+(1:nC),2] <- longDatlog$Y # continuous outcome
yB = yy[,1]
yC = yy[,2]

#######Add all survival covariates
# set up unique identifiers for the random-effects
longDatlog$idl <- ns+as.integer(longDatlog$id)
longDatlog$idl2 <- ns+ns+as.integer(longDatlog$id)
survDat$idl <- ns+as.integer(survDat$id)
survDat$idl2 <- ns+ns+as.integer(survDat$id)

cox_TRTs = as.factor(c(ifelse(survDat$trt=="0", "ref", "trt")))
cox_IntS = rep(1, length(cox_TRTs))
surv_inla_obj = inla.surv(time=survDat$deathTimes,event = survDat$d)
cox_ext = inla.coxph(surv_inla_obj ~  1+TRTs,
                     control.hazard=list(model="rw2", 
                     scale.model=TRUE,
                     diagonal=1e-4,
                     constr=F,
                     hyper=list(prec=list(prior="pc.prec",
                                          param=c(1,0.01)))),
                     data = c(list(surv_inla_obj = surv_inla_obj,
                                   TRTs = cox_TRTs,
                                   IDs = survDat$idl,
                                   IDsb = as.integer(survDat$id),
                                   IDs2 = survDat$idl2), as.list(survDat)))	
ns_cox = dim(cox_ext$data)[1] # for extended dataframe for poisson regression

###For other parts without survival part
# fixed effects
linear.covariate <- data.frame(
  InteB = c(rep(1,nB), rep(0,nC)), # intercept (binary part)
  InteC = c(rep(0,nB), rep(1,nC)), # intercept (continuous part)
  TIME = c(rep(0,nB),longDatlog$timej), # time (continuous part)
  TIMEb = c(longDat$timej,rep(0,nC)), # time (binary part)
  TRTc = c(rep(0,nB),as.numeric(longDatlog$trtY)-1), # treatment (continuous)
  TRTb = c(as.numeric(longDat$trtY)-1,rep(0,nC))) # treatment (binary)
# random-effects
random.covariate<-list(IDl=c(rep(NA,nB),longDatlog$idl), # random intercept (continuous)
                       IDb=c(as.integer(longDat$id),rep(NA,nC)), # random intercept (binary)
                       IDl2=c(rep(NA,nB),as.integer(longDatlog$idl2)), # random slope (continuous)
                       slopeCont=c(rep(NA,nB),longDatlog$timej)) # weight for random slope (continuous)
jointdf = data.frame(linear.covariate, random.covariate, yB, yC)
joint.data_cox <- c(as.list(inla.rbind.data.frames(jointdf, cox_ext$data)), 
                    cox_ext$data.list)
Yjoint = cbind(joint.data_cox$yB, joint.data_cox$yC, joint.data_cox$y..coxph) # outcomes
joint.data_cox$Y <- Yjoint

# conditional two-part joint model formula - update from the cox expansion
formulaJ= update(cox_ext$formula, Yjoint ~ . + InteB+InteC + TIME*TRTc+TIMEb*TRTb+
                   f(IDb, model="iid3d", n=3*ns,constr=F)+
                   f(IDl,copy="IDb")+
                   f(IDl2, slopeCont,copy="IDb")+
                   f(IDsb,copy="IDb",fixed=F)+
                   f(IDs,copy="IDl",fixed=F)+
                   f(IDs2,copy="IDl2",fixed=F))

#Fit model with INLA ()
TPinla <- inla(formulaJ,family = c("binomial", "gamma", cox_ext$family), 
               data=joint.data_cox,
               Ntrials=c(rep(1,length(longDat$Y)),rep(NA,nC),rep(NA,ns_cox)),
               control.predictor=list(compute=TRUE,link=1),#error gaussian
               E = joint.data_cox$E..coxph,
               control.family=list(list(control.link = list(model = "logit")),
                                   list(link='log',hyper = list(prec = list(initial = 2, fixed=FALSE))),
                                   list()),#variant = 1
               control.inla = list(strategy="adaptive"),
               control.fixed=list(remove.names="(Intercept)"),
               verbose=F)
print(summary(TPinla))
  