sink("s-CJS-Augmented.txt")
cat( "model { #model
# Priors
phi ~ dunif(0,1) # Survival (constant)
sigP ~ dunif(0,10) # Intercept in sigma estimate
sigP2<-sigP*sigP #variance of movement parameter sigma
lam0 ~ dgamma(0.1, 0.1) # Encounter rate
sigS~dunif(0,10)
tauXY <- 1/(sigS * sigS) #Variance for centroid spacing and placement

for (i in 1:C){ #m for all known individuals AT FIRST ENTRY
    SX[i,first[i]] ~ dunif(xl, xu)
    SY[i,first[i]] ~ dunif(yl, yu)
    z[i,first[i]] <- 1 # Known to be alive at entry into study
    

  for(j in 1:J){ #for all traps AT FIRST ENTRY  
    D2[i,j,first[i]] <- pow(SX[i,first[i]]-trapmat[j,1], 2) + pow(SY[i,first[i]]-trapmat[j,2],2)
    g[i,j,first[i]] <- lam0*exp(-D2[i,j,first[i]]/(2*sigP2))
    pmean[i,j,first[i]] <- 1- exp(-g[i,j,first[i]]) 
    tmp[i,j,first[i]] <- z[i,first[i]]*pmean[i,j,first[i]]
    y[i,j,first[i]] ~ dbin(tmp[i,j,first[i]], K)
  }
for (t in (first[i]+1):T) { #t #for all other time after
 #state process
    SX[i,t]~dnorm(SX[i, t-1], tauXY)T(xl,xu) # set priors for the X and Y coordinates of each individual
    SY[i,t]~dnorm(SY[i, t-1], tauXY)T(yl,yu)
    phiUP[i,t]<- z[i,t-1]*phi #animals dead of last time step cannot survive to the next, mu1?
    z[i,t] ~ dbern(phiUP[i,t]) #probability of being alive 
  for(j in 1:J){
#observation process 
    D2[i,j,t] <- pow(SX[i,t]-trapmat[j,1], 2) + pow(SY[i,t]-trapmat[j,2],2)
    g[i,j,t] <- lam0*exp(-D2[i,j,t]/(2*sigP2)) #gaussian hazard detection 
    pmean[i,j,t] <- 1- exp(-g[i,j,t]) #probabiity of capture
    tmp[i,j,t] <- z[i,t]*pmean[i,j,t] #p.eff
    y[i,j,t] ~ dbin(tmp[i,j,t], K)
  }           
} # t
}#M

for (i in (C+1):M){ #m for all UNKNOWN individuals AT FIRST ENTRY
    z[i,first[i]] ~ dbern(phi) #Augmented individual not seen and estimated if alive
    SX[i,first[i]] ~ dunif(xl, xu)
    SY[i,first[i]] ~ dunif(yl, yu)

  for(j in 1:J){ #for all traps AT FIRST ENTRY  
    D2[i,j,first[i]] <- pow(SX[i,first[i]]-trapmat[j,1], 2) + pow(SY[i,first[i]]-trapmat[j,2],2)
    g[i,j,first[i]] <- lam0*exp(-D2[i,j,first[i]]/(2*sigP2))
    pmean[i,j,first[i]] <- 1- exp(-g[i,j,first[i]]) 
    tmp[i,j,first[i]] <- z[i,first[i]]*pmean[i,j,first[i]]
    y[i,j,first[i]] ~ dbin(tmp[i,j,first[i]], K)
  }
for (t in (first[i]+1):T) { #t #for all other time after
 #state process
    SX[i,t]~dnorm(SX[i, t-1], tauXY)T(xl,xu) # set priors for the X and Y coordinates of each individual
    SY[i,t]~dnorm(SY[i, t-1], tauXY)T(yl,yu)
    phiUP[i,t]<-phi #*z[i,t-1] #animals dead of last time step cannot survive to the next, mu1? But maybe should be just phi?
     z[i,t] ~ dbern(phiUP[i,t]) #probability of being alive 
  for(j in 1:J){
#observation process 
    D2[i,j,t] <- pow(SX[i,t]-trapmat[j,1], 2) + pow(SY[i,t]-trapmat[j,2],2)
    g[i,j,t] <- lam0*exp(-D2[i,j,t]/(2*sigP2)) #gaussian hazard detection 
    pmean[i,j,t] <- 1- exp(-g[i,j,t]) #probabiity of capture
    tmp[i,j,t] <- z[i,t]*pmean[i,j,t] #p.eff
    y[i,j,t] ~ dbin(tmp[i,j,t], K)
  }           
} # t
}#M

N <- sum(z[1:M,T])
D <- N/144  #(xlim[2]*ylim[2])
} #model

 ",fill=TRUE)  
sink()


#Set up to run model-----------------------------------------------------------     

#set up intial values for z

zin<-matrix(1, dim(simdat$y)[1], T)
for(i in 1:dim(simdat$y)[1]){
  for(t in 1:simdat$first[i]){
    zin[i,t]<-NA
  }}

zin_aug <- matrix(NA, (M-dim(simdat$y)[1]), T)

zin <- rbind(zin, zin_aug)

#Augment observational array
library(abind)
y = simdat$y

augmented_row <- array(0, dim = c((M-dim(simdat$y)[1]), dim(simdat$y)[2], dim(simdat$y)[3])) #blank observation rows for M-N individuals
augmented_row
# Combine using abind
y<- abind(y, augmented_row, along = 1)

#Augment First Seen Vector

unseen <- rep(1, (M-dim(simdat$y)[1]))

first <- c(simdat$first, unseen)

#function to creat initial S for each year
Sin<-function(first=first,T=T, M=M, xl=xl, xu=xu, yl=yl, yu=yu,ntot=ntot){
  SX<-SY<-matrix(NA, nrow=ntot, ncol=T)
  for(i in 1:ntot){
    for (t in first[i]:T){
      SX[i,t]<-runif(1, xl, xu)
      SY[i,t]<-runif(1, yl, yu)
      traps<-which(y[i,,t]>0)
      if(length(traps)>0){
        SX[i,t]<- mean(grid[traps,1])
        SY[i,t]<- mean(grid[traps,2])
      }
    }
  }
  return(list(SX, SY))
}
test <- Sin(first=first,T=T, xl=xl, xu=xu, yl=yl, yu=yu, ntot=dim(y)[1])[[1]]
test

#statespace is 4x, as is data generation
data.leah<-list(y=y, first=first,T=5,xl=xl, xu=xu, yl=yl, yu=yu,K=K,
                M=dim(y)[1],J=J, trapmat=as.matrix(grid), C = as.numeric(dim(simdat$y)[1]))

inits.leah = function() {list(phi=runif(1, .5, 1),sigP=runif(1,1,2), lam0=runif(1), z=zin, sigS=runif(1,0.1,3),
                              SX=Sin(first=first,T=T, xl=xl, xu=xu, yl=yl, yu=yu, ntot=dim(y)[1])[[1]],
                              SY=Sin(first=first,T=T, xl=xl, xu=xu, yl=yl, yu=yu, ntot=dim(y)[1])[[2]]
) }
params<-c("phi","sigP", "sigS", "lam0", "N", "D", "z")
#cjs.mod2 <- textConnection(cjs.mod)
#mod<-jags.model(cjs.mod, data, inits, n.chains=3, n.adapt=100)
mod.out.aug <- jags.model("s-CJS-Augmented.txt", data.leah, inits.leah, n.chains=3)
out.aug<-coda.samples(mod.out.aug, params, n.iter=10000, thin=3)
mod.out.aug$

library(MCMCvis)
library(mcmcplots)    
MCMCplot(out.aug, params = c('phi', 'N', 'lam0', 'sigP', 'sigS'))

mcmcplots::mcmcplot(out.aug)

LINE.out <- coda.samples(out.aug, c("z"), n.iter=1000)
out.aug$mcmc.list
as.matrix(out.aug[1])

summary(out.aug)
