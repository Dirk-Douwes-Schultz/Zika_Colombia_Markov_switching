#code to run the zero-state space-time Markov switching model from the paper
#will run the model on the Zika data and check convergence
rm(list=ls())

# Packages ----------------------------------------------------------------
library(nimble)
library(tidyverse)
library(dlnm)
library(splines)
library(coda)
library(parallel) 
library(doParallel) 
# memory.limit(size=1E10)
# 
setwd("~/GitHub/Zika_Colombia_Markov_switching")


# Data --------------------------------------------------------------------

load("zika_cases_covariates_colombia.RData") 
all <- zika
rm(zika)

#cases_it=y_it
zika <- all %>% 
  filter(epiweek >= 22 & epiyear == 2015 | #1st reported cases in week 26/2015
           epiyear == 2016) %>% 
  replace_na(list(zika = 0)) 

nareas <- length(unique(zika$cod_mun_res))
nweeks <- length(unique(paste0(zika$epiweek,zika$epiyear)))
cases <- matrix(zika$zika, nrow = nareas, ncol = nweeks, byrow = TRUE)

#cumulative incidence_i t-4
cuminc <- all %>% 
  filter(epiweek >=18 & epiyear == 2015 | epiweek <= 35 & epiyear == 2016) %>% 
  replace_na(list(zika = 0)) %>% 
  group_by(cod_mun_res) %>% 
  mutate(cumcases = cumsum(zika),
         cuminc = (cumcases/pop)*10000) 
cuminc <- matrix(cuminc$cuminc, nrow = nareas, ncol = nweeks, byrow = TRUE)

#Rain_i t-4
Rain <- all %>% 
  filter(epiweek >=18 & epiyear == 2015 | epiweek <= 35 & epiyear == 2016)
Rain <- matrix(Rain$precip, nrow = nareas, ncol = nweeks, byrow = TRUE)
Rain <- (Rain - mean(Rain))/sd(Rain)

#Humidity_i t-1
hum <- all %>% 
  filter(epiweek >=21 & epiyear == 2015 | epiweek <= 38 & epiyear == 2016)
hum <- matrix(hum$relhum, nrow = nareas, ncol = nweeks, byrow = TRUE)
hum <- (hum - mean(hum))/sd(hum)

#Temp_i t-1
Temp <- all %>% 
  filter(epiweek >=21 & epiyear == 2015 | epiweek <= 38 & epiyear == 2016)
#filter(epiweek >= 26 & epiyear == 2015 | epiyear == 2016)
Temp <- matrix(Temp$tmax, nrow = nareas, ncol = nweeks, byrow = TRUE)
Temp <- (Temp - mean(Temp))/sd(Temp)

#NDVI_i 
covsfixed <- all %>% 
  filter(epiweek >= 22 & epiyear == 2015 | epiyear == 2016) %>% 
  group_by(cod_mun_res) %>% 
  summarise(across(c(ndviterra, pop, popdens_km2), mean, na.rm = TRUE),
            across(c(elevation_m,nbi2018), unique))
ndvi <- covsfixed$ndviterra
ndvi <- (ndvi - mean(ndvi))/sd(ndvi)

#N_i 
N <- covsfixed$pop

#NBI_i
nbi <- covsfixed$nbi2018 
nbi <- (nbi - mean(nbi))/sd(nbi)

# population density_i
pop <- log(covsfixed$popdens_km2)
pop <- (pop - mean(pop))/sd(pop)

#elevation_i
elevation <- covsfixed$elevation_m 
elevation <- (elevation - mean(elevation))/sd(elevation)

# Department
dpt <- all %>% 
  select(DPTO,NOMBRE_DPT,cod_mun_res) %>% 
  unique() 

dpt$dpto <- factor(dpt$DPTO, labels = seq(1:33))
dpto <- as.numeric(dpt$dpto)

#also need S_it=NA if y_it=0 and S_it=1 if y_it>0
S <- matrix(nrow=nareas,ncol=nweeks)
for(i in 1:nareas){
  for(t in 1:nweeks){
    if(cases[i,t]>0){
      S[i,t] <- 3
    }
  }
}


# Model -------------------------------------------------------------------

dengeeConsts <- list(N=nareas,
                     T=nweeks, 
                     pop=pop,
                     Temp=Temp,
                     Rain=Rain,
                     dpto = dpto,
                     ndvi=ndvi,
                     elevation=elevation,
                     nbi=nbi,
                     hum=hum,
                     cuminc=cuminc,
                     lcuminc = log(cuminc+1),
                     mlcuminc = mean(log(cuminc+1)),
                     count=count,
                     psi=cases,
                     mpsi=mean(cases),
                     lpsi=log(cases+1),
                     mlpsi=mean(log(cases+1)))

dengeeData <- list(y=cases,S=S) 

dengeeCode <- nimbleCode({
  
  #priors
  for (i in 1:33) {
    beta0[i] ~ dnorm(beta0cent, precision_beta0)
    alpha1[i] ~ dnorm(alphacent1, precision_alpha1)
    alpha6[i] ~ dnorm(alphacent6, precision_alpha6)
    alpha5[i] ~ dnorm(alphacent5, precision_alpha5)
  }
  beta0cent ~ dnorm(0, sd = 5)
  alphacent1 ~ dnorm(0, sd = 5)
  alphacent6 ~ dnorm(0, sd = 5)
  alphacent5 ~ dnorm(0, sd = 5)
  precision_beta0 ~ dgamma(.1,.1)
  precision_alpha1 ~ dgamma(.1,.1)
  precision_alpha6 ~ dgamma(.1,.1)
  precision_alpha5 ~ dgamma(.1,.1)
  
  #beta0 ~ dnorm(0, sd=100)
  beta1 ~ dnorm(0, sd=100)
  beta2 ~ dnorm(0, sd=100)
  beta3 ~ dnorm(0, sd=100)
  beta4 ~ dnorm(0, sd=100)
  beta5 ~ dnorm(0, sd=100)
  beta6 ~ dnorm(0, sd=100)
  beta7 ~ dnorm(0, sd=100)
  beta8 ~ dnorm(0, sd=100)
  beta9 ~ dnorm(0, sd=100)
  
  rho ~ dnorm(0, sd=100)
  r2 ~ dunif(0,50)
  
  for(j in 2:4){
    alpha[j] ~ dnorm(0, sd=100)
  }
  for(j in 8:12){
    alpha[j] ~ dnorm(0, sd=100)
  }
  #alpha[6] ~ dnorm(0, sd = 5)
  alpha[7] ~ dnorm(0, sd = 1.82)
  
  for(i in 1:N){
    b0[i] ~ dnorm(beta0[dpto[i]]+
                    beta3*pop[i] +
                    beta4*ndvi[i] +
                    beta5*elevation[i]+
                    beta6*nbi[i],  
                  precision_b0)
    b[i] ~ dnorm(rho, precision_b)
  }
  precision_b0 ~ dgamma(.1,.1)
  precision_b ~ dgamma(.1,.1)
  
  #likelihood
  indi[1:3] <- c(0,0,1)
  for(i in 1:N) {
    for(t in 2:T){
      mup[i,t] <- exp(b0[i] + 
                        beta1*Rain[i,t] +
                        beta2*Temp[i,t] +
                        beta7*hum[i,t] +
                        beta8*(lcuminc[i,t]-mlcuminc) + 
                        beta9*((lcuminc[i,t]-mlcuminc)^2)) *
        psi[i,t-1]+exp(b[i])
      
      p[i,t] <- r2/(r2+(indi[S[i,t]])*mup[i,t]) - 1e-10*(1-indi[S[i,t]])
      y[i,t] ~ dnegbin(p[i,t],r2)
      #mu[i,t] <- mup[i,t]*(S[i,t])+.00001
      #y[i,t] ~ dpois(mu[i,t])
      # like[i,t] <- dpois(x=y[i,t],lambda=mu[i,t])
      # yfit[i,t] ~ dpois(mu[i,t])
    }
  }
  #markov chain
  init_probs[1:3] <- c(.95,0,.05)
  for(i in 1:N) {
    S[i,1] ~  dcat(init_probs[1:3])
    for(t in 2:T){ 
      
      #calc transition matrix
      tm[i,t,1,1:3] <- c(1-p13[i,t],0,p13[i,t])
      tm[i,t,2,1:3] <- c(0,1-p23[i,t],p23[i,t])
      tm[i,t,3,1:3] <- c(0,1-p33[i,t],p33[i,t])
      
      logit(p13[i,t]) <- alpha1[dpto[i]]+
        alpha[2]*pop[i] +
        alpha[3]*Temp[i,t] +
        alpha[4]*Rain[i,t] +
        alpha[8]*ndvi[i] +
        alpha[9]*elevation[i] +
        alpha[10]*nbi[i] +
        alpha[11]*hum[i,t] 
      
      logit(p23[i,t]) <- alpha5[dpto[i]]+
        alpha[2]*pop[i] +
        alpha[3]*Temp[i,t] +
        alpha[4]*Rain[i,t] +
        alpha[8]*ndvi[i] +
        alpha[9]*elevation[i] +
        alpha[10]*nbi[i] +
        alpha[11]*hum[i,t]
      
      logit(p33[i,t]) <- alpha6[dpto[i]]+
        alpha[7]*(lpsi[i,t-1]-mlpsi) +
        alpha[12]*pop[i]
      
      S[i,t] ~ dcat(tm[i,t,S[i,t-1],1:3])
    }
  }
})

inits <- list(beta0 = rnorm(n=33,mean=0,sd=.1),
              beta0cent = rnorm(n=1,mean=0,sd=.01),
              beta1 = rnorm(n=1,mean=0,sd=.01),
              beta2 = rnorm(n=1,mean=0,sd=.01),
              beta3 = rnorm(n=1,mean=0,sd=.01),
              beta4 = rnorm(n=1,mean=0,sd=.01),
              beta5 = rnorm(n=1,mean=0,sd=.01),
              beta6 = rnorm(n=1,mean=0,sd=.01),
              beta7 = rnorm(n=1,mean=0,sd=.01),
              beta8 = rnorm(n=1,mean=0,sd=.01),
              beta9 = rnorm(n=1,mean=0,sd=.01),
              precision_b0 = 1, 
              precision_b = 1,
              precision_beta0 = 1, 
              precision_alpha1 = 1,
              precision_alpha6 = 1, 
              precision_alpha5 = 1, 
              rho = rnorm(n=1,mean=0,sd=1) ,
              S = matrix(rep(3,nareas*nweeks),nrow = nareas, ncol = nweeks),
              alpha = rnorm(n=12,mean=0,sd=c(1,.01,1,.1,1,1,1,.1,.1,.1,.1,.01)),
              alpha1 = rnorm(n=33,mean=0,sd=1),
              alpha6 = rnorm(n=33,mean=0,sd=1),
              alpha5 = rnorm(n=33,mean=0,sd=1),
              alphacent1 = rnorm(n=1,mean=0,sd=.01),
              alphacent6 = rnorm(n=1,mean=0,sd=.01),
              alphacent5 = rnorm(n=1,mean=0,sd=.01),
              b = rnorm(n=nareas,mean=0,sd=.01),
              b0 = rnorm(n=nareas,mean=0,sd=.01),
              r2=runif(n=1,min=0,max=20))

dengeemodel <- nimbleModel(dengeeCode, dengeeConsts, dengeeData, inits)

Cdengee <- compileNimble(dengeemodel)


#now it drafts and tests a FFBS sampler
#testing
#model$getLogProb()

#arguments
model <- dengeemodel
dq1 <- !dengeemodel$isData("S[1, ]")
target <- "S[1, dq1]"

model$S[1, dq1]

#setup
nnames <- model$expandNodeNames(target)
times <- as.numeric(gsub(".*\\[.*,(.*)\\].*", "\\1", nnames))
loc <- as.numeric(gsub(".*\\[(.*),.*\\].*", "\\1", nnames[[1]]))
numnodes <- length(nnames)
startbs <- max(times)
startfs <- min(times)
Mt <- length(model$y[loc, ])
startmidfzt <- ifelse(startfs==1,2,1)
endmidfzt <- ifelse(startbs==Mt,numnodes-1,numnodes)
startizt <- ifelse(startbs==Mt,2,1)
#all dependencies together, not same as "dependencies" because that deletes duplicates
calcNodes <- model$getDependencies(target)
#grab the dependencies of each target node
dependencies <- NULL
start_depend <- rep(NA,numnodes)
start_depend[1] <- 1
end_depend <- rep(NA,numnodes)
index <- 1
for (n in nnames){
  d <- model$getDependencies(n)
  dependencies <- c(dependencies,d)
  end_depend[index] <- length(d)+start_depend[index]-1
  start_depend[index+1] <- end_depend[index]+1
  index <- index+1
}
#store filtered probabilities
q <- matrix(nrow=numnodes,ncol=3)
q[1,1:3] <- c(-.99,-.99,-.99)
#log likelihood
ll <- matrix(nrow=numnodes,ncol=3)
ll[1,1:3] <- c(-.99,-.99,-.99)

#now run 
#start with ct=1
#in a non coupled filter it is just the initial state distribution
q[1,1:3] <- model$init_probs

for(zt in 2:numnodes){
  
  ct <- times[zt]
  #check if previous time is data or not
  if(times[zt-1]==ct-1){
    q_tm1 <- q[zt-1,1:3]
  }else{
    q_tm1 <- c(0,0,1)
  }
  
  #predictive probability
  p <- t(model$tm[loc,ct,1:3,1:3]) %*% asCol(q_tm1[1:3])
  
  #log likelihoods, remember y[loc,ct]=0
  ll[zt,1] <- 0
  ll[zt,2] <- 0
  ll[zt,3] <- dnbinom(x=0,size = model$r2[1] ,
                      prob = model$r2[1]/(model$r2[1]+model$mup[loc,ct]),log =TRUE)
  
  nl <- ll[zt,1:3]+log(p[1:3,1])
  nls <- nl-max(nl)
  q[zt,1:3] <- exp(nls)/sum(exp(nls))
  
}

#testing
# model$S[1, dq1]
# model$getLogProb()
# model$calculate()

#now backwards sampling

#start with Mt 
if(startbs==Mt){
  prev <- model$S[loc,Mt]
  model$S[loc,Mt] <- rcat(1,q[numnodes,1:3])
  if(model$S[loc,Mt]!=prev){
    model$calculate(nodes = dependencies[start_depend[numnodes]:end_depend[numnodes]])
  }
}

for(izt in startizt:(numnodes)){
  zt <- numnodes-izt+1
  ct <- times[zt]
  prev <- model$S[loc,ct]
  fs <- model$S[loc,ct+1]
  trans <- model$tm[loc,ct+1,1:3,fs]
  lp <- log(trans)+log(q[zt,1:3])
  lp <- lp-max(lp)
  b <- exp(lp)/sum(exp(lp))
  news <- rcat(1,b)
  model$S[loc,ct] <- news
  if(model$S[loc,ct]!= prev){
    model$calculate(nodes = dependencies[start_depend[zt]:end_depend[zt]])
  }
}

#testing
# model$S[1, dq1]
# model$getLogProb()
# model$calculate()

FFBS <- nimbleFunction(
  
  contains = sampler_BASE,
  
  setup = function(model, mvSaved, target, control) {
    
    #setup
    nnames <- model$expandNodeNames(target)
    times <- as.numeric(gsub(".*\\[.*,(.*)\\].*", "\\1", nnames))
    loc <- as.numeric(gsub(".*\\[(.*),.*\\].*", "\\1", nnames[[1]]))
    numnodes <- length(nnames)
    startbs <- max(times)
    startfs <- min(times)
    Mt <- length(model$y[loc, ])
    startmidfzt <- ifelse(startfs==1,2,1)
    endmidfzt <- ifelse(startbs==Mt,numnodes-1,numnodes)
    startizt <- ifelse(startbs==Mt,2,1)
    #all dependencies together, not same as "dependencies" because that deletes duplicates
    calcNodes <- model$getDependencies(target)
    #grab the dependencies of each target node
    dependencies <- NULL
    start_depend <- rep(NA,numnodes)
    start_depend[1] <- 1
    end_depend <- rep(NA,numnodes)
    index <- 1
    for (n in nnames){
      d <- model$getDependencies(n)
      dependencies <- c(dependencies,d)
      end_depend[index] <- length(d)+start_depend[index]-1
      start_depend[index+1] <- end_depend[index]+1
      index <- index+1
    }
    #store filtered probabilities
    q <- matrix(nrow=numnodes,ncol=3)
    q[1,1:3] <- c(-.99,-.99,-.99)
    #log likelihood
    ll <- matrix(nrow=numnodes,ncol=3)
    ll[1,1:3] <- c(-.99,-.99,-.99)
    
    
    
  },
  
  
  run = function() {
    
    #now run 
    #start with ct=1
    #in a non coupled filter it is just the initial state distribution
    q[1,1:3] <<- model$init_probs
    
    for(zt in 2:numnodes){
      
      ct <- times[zt]
      #check if previous time is data or not
      if(times[zt-1]==ct-1){
        q_tm1 <- q[zt-1,1:3]
      }else{
        q_tm1 <- c(0,0,1)
      }
      
      #predictive probability
      p <- t(model$tm[loc,ct,1:3,1:3]) %*% asCol(q_tm1[1:3])
      
      #log likelihoods, remember y[loc,ct]=0
      ll[zt,1] <<- 0
      ll[zt,2] <<- 0
      ll[zt,3] <<- dnbinom(x=0,size = model$r2[1] ,
                           prob = model$r2[1]/(model$r2[1]+model$mup[loc,ct]),log =TRUE)
      
      nl <- ll[zt,1:3]+log(p[1:3,1])
      nls <- nl-max(nl)
      q[zt,1:3] <<- exp(nls)/sum(exp(nls))
      
    }
    
    #now backwards sampling
    
    #start with Mt 
    if(startbs==Mt){
      prev <- model$S[loc,Mt]
      model$S[loc,Mt] <<- rcat(1,q[numnodes,1:3])
      if(model$S[loc,Mt]!=prev){
        model$calculate(nodes = dependencies[start_depend[numnodes]:end_depend[numnodes]])
      }
    }
    
    for(izt in startizt:(numnodes)){
      zt <- numnodes-izt+1
      ct <- times[zt]
      prev <- model$S[loc,ct]
      fs <- model$S[loc,ct+1]
      trans <- model$tm[loc,ct+1,1:3,fs]
      lp <- log(trans)+log(q[zt,1:3])
      lp <- lp-max(lp)
      b <- exp(lp)/sum(exp(lp))
      news <- rcat(1,b)
      model$S[loc,ct] <<- news
      if(model$S[loc,ct]!= prev){
        model$calculate(nodes = dependencies[start_depend[zt]:end_depend[zt]])
      }
    }
    
    copy(from = model, to = mvSaved, row = 1, 
         nodes = calcNodes, logProb = TRUE)
    
  },
  
  methods = list(   reset = function () {}   )
  
)

dengeeConf <- configureMCMC(dengeemodel, print = TRUE)


#add FFBS sampler, this might not work in new versions of nimble
for(loc in 1:1121){
  
  
  dq <- !dengeemodel$isData(paste0("S[",loc,", ]"))
  if(sum(dq)>0){
    dengeeConf$removeSampler(paste0("S[",loc,", dq]"))
    dengeeConf$addSampler(target = dengeemodel$expandNodeNames(paste0("S[",loc,", dq]")),
                          type = "FFBS")
  }
  
}

print(dengeeConf)

dengeeConf$addMonitors(c("S","b","b0","beta0","alpha1","alpha6","alpha5"))


dengeeMCMC <- buildMCMC(dengeeConf)

CdengeeMCMC <- compileNimble(dengeeMCMC, project = dengeemodel, resetFunctions = TRUE)


initsFunction <- function() list(beta0 = rnorm(n=33,mean=0,sd=.1),
                                 beta0cent = rnorm(n=1,mean=0,sd=.01),
                                 beta1 = rnorm(n=1,mean=0,sd=.01),
                                 beta2 = rnorm(n=1,mean=0,sd=.01),
                                 beta3 = rnorm(n=1,mean=0,sd=.01),
                                 beta4 = rnorm(n=1,mean=0,sd=.01),
                                 beta5 = rnorm(n=1,mean=0,sd=.01),
                                 beta6 = rnorm(n=1,mean=0,sd=.01),
                                 beta7 = rnorm(n=1,mean=0,sd=.01),
                                 beta8 = rnorm(n=1,mean=0,sd=.01),
                                 beta9 = rnorm(n=1,mean=0,sd=.01),
                                 precision_b0 = 1, 
                                 precision_b = 1,
                                 precision_beta0 = 1, 
                                 precision_alpha1 = 1,
                                 precision_alpha6 = 1, 
                                 precision_alpha5 = 1, 
                                 rho = rnorm(n=1,mean=0,sd=1) ,
                                 S = matrix(rep(3,nareas*nweeks),nrow = nareas, ncol = nweeks),
                                 alpha = rnorm(n=12,mean=0,sd=c(1,.01,1,.1,1,1,1,.1,.1,.1,.1,.01)),
                                 alpha1 = rnorm(n=33,mean=0,sd=1),
                                 alpha6 = rnorm(n=33,mean=0,sd=1),
                                 alpha5 = rnorm(n=33,mean=0,sd=1),
                                 alphacent1 = rnorm(n=1,mean=0,sd=.01),
                                 alphacent6 = rnorm(n=1,mean=0,sd=.01),
                                 alphacent5 = rnorm(n=1,mean=0,sd=.01),
                                 b = rnorm(n=nareas,mean=0,sd=.01),
                                 b0 = rnorm(n=nareas,mean=0,sd=.01),
                                 r2=runif(n=1,min=0,max=20))


samples <- runMCMC(CdengeeMCMC, 
                   niter = 400000,
                   nchains = 3, 
                   nburnin = 100000 ,
                   samplesAsCodaMCMC = TRUE, 
                   thin = 75, 
                   inits = initsFunction)


parameter_vector <- c(paste0("beta0[",1:33,"]"),paste0("alpha1[",1:33,"]"),
                      paste0("alpha5[",1:33,"]"),paste0("alpha6[",1:33,"]"),
                      paste0("alpha[",c(2,3,4,7,8,9,10,11,12),"]"),
                      paste0("beta",2:9),
                      "beta0cent","alphacent1","alphacent6","alphacent5",
                      "precision_beta0","precision_alpha1","precision_alpha6","precision_alpha5",
                      "precision_b","precision_b0","rho")

#check convergence
gelman.diag(samples[,parameter_vector])
min(effectiveSize(samples[,parameter_vector]))

#test the municipal intercepts
AR_gr <- rep(NA,1121)
BL_gr <- rep(NA,1121)

for(i in 1:1121){
  AR_gr[i] <- gelman.diag(samples[,paste0("b0[",i,"]")])$psrf[1,2]
  BL_gr[i] <- gelman.diag(samples[,paste0("b[",i,"]")])$psrf[1,2]
}

max(AR_gr)
max(BL_gr)

#trace plots
for(p in 1:53){
  j <- (3*p)-2
  plot(samples[,parameter_vector[j:(j+2)]])
}


plot(samples[,parameter_vector[159:160]])


