library(jagsUI) 
library(MCMCvis)
library(R2jags)



#Simulation -------------------------------------------------------------------


#Define Parameters-----------------------------------

months <- 3 #number of months in dataset TOTAL
primary_occ <- c(1:3) #number actually trapped, skip # if not trapped
secondary_occ <- c(rep(3,3)) # number of secondary occasions PER MONTH, list 0 if not trapped

ind <- rep(20, months-1) #number of individuals marked per month, list format

#currently written to assume this will be consistent, but it can be written to vary
phi <- rep(0.65, months-1) #survival
p <- rep(0.3, months) #capture probability
c <- rep(0.4, months) #secondary capture probability, assumes more recapture upon initial capture

#Make Matrices---------------------------------------
PHI <- matrix(phi, ncol = months-1, nrow = sum(ind))
P <- matrix(p, ncol=months, nrow=sum(ind))
C <- matrix(c, ncol=months, nrow=sum(ind))

#define function to simulate a capture history matrix (CH)
simul.cjs.rb <- function(PHI, P, C, marked, n.months, prim.occasions, n.sec.occasions){
  z <- array(0, dim = c(sum(ind), months)) # z is actual state
  y <- array(0, dim = c(sum(ind), max(secondary_occ), months)) # y is Ch observed so includes secondary occasions
  
  
  #define a vector with the occasion of marking
  marked.occ <- rep(1:length(ind), ind[1:length(ind)])
  
  #fill the CH Matrix
  for(i in 1:sum(ind)){
    z[i, marked.occ[i]] <- 1		#put a 1 at the release occasion
    
    # for first month caught 
    #Bernouli trial: is indiv captured?
    ########  secondary occasions, d for days
    for(d in 1:max(secondary_occ)){
      p.eff <- ifelse(sum(y[i, 1:(d-1), marked.occ[i]])==0, P[i, marked.occ[i]], C[i, marked.occ[i]]) #if caught any time previously in this session then use c instead of p
      y[i, d, marked.occ[i]] <- rbinom(1, 1, prob = p.eff)
    } #d
    # if never caught then randomly pick a secondary occasion for capture
    if(sum(y[i, , marked.occ[i]],na.rm=TRUE)==0){y[i, sample(1:secondary_occ[marked.occ[i]], 1), marked.occ[i]] <- 1}
    
    if(marked.occ[i] == months) next	#starts next iter of loop if caught only at last occasion
    
    for(m in (marked.occ[i]+1):months){ # m is primary occasion (month)
      #p.eff <- array(NA, dim = sum(marked))
      mu1 <- PHI[i, m-1] * z[i, m-1] # this assures that animals stay dead
      z[i, m] <-  rbinom(1, 1, mu1) 		#Bernouli trial: does animal survive
      
      if(mu1==0) break				# if dead, move to next indiv
      
      #Bernouli trial: is indiv captured?
      ########  secondary occasions, d for days
      for(d in 1:max(secondary_occ)){
        p.eff <- ifelse(sum(y[i, 1:(d-1), m])==0, P[i, m], C[i, m]) * z[i, m] #if caught any time previously in this session (m) then c, and if not alive, can't be caught
        y[i, d, m] <- rbinom(1, 1, prob = p.eff)
      } #d
    } #m
  } #i
  
  ## put NAs in for days not trapped  
  inds <- which(n.sec.occasions<max(n.sec.occasions))
  for(m in inds){
    y[ , (n.sec.occasions[m]+1):max(n.sec.occasions) , m] <- NA
  }
  ## put NAs in for months not trapped  
  inds <- which(1 - (1:n.months %in% prim.occasions )==1)
  for(m in inds){
    y[ , , m] <- NA
  }
  
  return(list(true.state=z,observed=y))	
}

sim.data=simul.cjs.rb(PHI, P, C, marked, months, primary_occ, secondary_occ)

sim.data

y.observed<- sim.data$observed
# remove indivdiuals never caught (because of incomplete sampling)
y.observed <- y.observed[(which(apply(y.observed,1,sum, na.rm = TRUE)!=0)) , , ]
y.observed
#Organize Data------------------------------------------------------------------


# Remove animals only caught at last primary occasion (don't help and mess with code)

CH.primary <- apply(y.observed,c(1,3),sum,na.rm=TRUE)
CH.primary[,which(1 - (1:months %in% primary_occ )==1)] <- NA # make sure keep NAs when whole month not trapped 
CH.primary <- replace(CH.primary,CH.primary>1,1)

CH.primary

# create a vector of first marking
get.first <- function(x) min(which(x!=0))
f <- apply(CH.primary, 1, get.first)
f
#delete individuals that were caught for the first time at the last primary occasion
# they don't provide data and mess with code
rms <- which(f==dim(CH.primary)[2])
if(length(rms)>0){
  CH.secondary <- CH.secondary[-rms, , ]
  CH.primary <- CH.primary[-rms, ]
  f <- f[-rms]
}
f

# create a separate f in terms of months captured (not full months)
Prim <- match(1:months, primary_occ)

f.prim <- Prim[f]

f.prim

f

#function to create matrix with info about known latent state z
known.state.cjs=function(ch){
  state=ch
  for(i in 1:dim(ch)[1]){
    n1=min(which(ch[i,]==1))
    n2=max(which(ch[i,]==1))
    state[i,n1:n2]=1
    state[i,n1]=NA			#only filling in those that were 0s but we know were alive because caught before and after
  }
  state[state==0]=NA
  return(state)
}

known.state.cjs(CH.primary)

#Prime Data for Model---------------------------------------------------------------

# [indiv, secondary occasions, primary occasions]
##### Bundle data
bugs.data=list(y = y.observed,
               f = f,
               f.prim = f.prim,
               nind = dim(y.observed)[1],
               n.secondary.occasions = secondary_occ, # contains 0's for months not trapped
               n.primary.occasions = months, # total span of months (incl not trapped)
               prim.occasions = primary_occ, # months out of n.months actually trapped (if all trapped then same as previous line)
               z = known.state.cjs(CH.primary))

###### function to create matrix of initial values for latent state z
# we shouldn't give initial values for those elements of z whose value is specified in the data.
# they get an NA
cjs.init.z <- function(ch,f,prim.occasions){
  ch <- replace(ch,is.na(ch),0)
  for(i in 1:dim(ch)[1]){
    if(sum(ch[i,],na.rm=TRUE)==1) next
    n2=max(which(ch[i,]==1))
    ch[i,f[i]:n2]=NA
  }
  for(i in 1:dim(ch)[1]){
    ch[i,1:f[i]]=NA
  }
  return(ch)	
}

cjs.init.z(CH.primary,f)



#initial values
inits=function(){list(
  z=cjs.init.z(CH.primary,f),
  mean.phi=runif(1,0,1),
  mean.p=runif(1,0,1),
  mean.c=runif(1,0,1))}

#parameters monitored
parameters=c("mean.phi","mean.p","mean.c")

#MCMCsettings
ni=10000
nt=6
nb=5000
nc=3

#CJS Robust Model----------------------------------------------------------------

sink("robust_cjs.bug")
cat("					######<--------------------- uncomment 
model{
	
###############Priors and constraints
mean.phi ~ dunif(0, 1) # prior for phi
mean.p ~ dunif(0, 1)   # prior for p
mean.c ~ dunif(0, 1)   # prior for c
    
for(i in 1:nind){
	
	for(m in f[i]:n.primary.occasions){  ### for p need every time for phi need -1
		# phi has only 2 dimensions [indiv, and primary occasions]
    phi[i,m] <- mean.phi   # could specify covariates here
    # p and c have 3 dimensions [indiv, secondary, primary]		
    
    for(d in 1:max(n.secondary.occasions)){
      p[i, d, m] <- mean.p  # could specify temporal or individual covariates here
      c[i, d, m] <- mean.c  # ditto
	   } #d for days
    } #m for months
	} #i for individual

#############Likelihood 		

for(i in 1:nind){
	#define latent state at first capture 
  # dimensions [individual, primary session (month)]
	z[i,f[i]] <- 1		# z is true (latent) state alive or dead, know alive at first capture
	
	for(m in  (f[i]+1):n.primary.occasions){  
		
	  #state process				# alive or dead
	  z[i, m] ~ dbern(mu1[i, m]) 		#mu1 is probability alive
		mu1[i, m] <- phi[i, m] * z[i, m-1] # this assures that animals stay dead
	} # m
	
	for(m in prim.occasions[f.prim[i]:length(prim.occasions)] ) { 
	  #observation process			# caught or not
	  #first secondary occasion within a primary occasion:
	  y[i, 1, m] ~ dbern(p.eff[i, 1, m])
	  p.eff[i, 1, m] <- z[i, m] * p[i, 1, m]   
	  
	  #loop over rest of secondary occasions per primary
	  for(d in 2:n.secondary.occasions[m]){ ## add index for number of secondary occasions if not consistent
		  y[i, d, m] ~ dbern(p.eff[i, d, m]) 		# p.eff is prob of capture
		  p.eff[i, d, m] <- z[i, m] * ifelse(sum(y[i, 1:(d-1), m])==0, p[i, d, m], c[i, d, m])	
      # capture prob= (p if not caught previously that session or c if was caught that session) 
      # multiply by if it was alive (so can't capture animal that's not alive)
		 
	   } #d
		} #m
	} #i
}
",fill=TRUE)  #####<----------------uncomment this
sink()



#Run Model----------------------------------------------------------------------
library(jagsUI)

date()
## Call JAGS from R
robust.cjs=jags(bugs.data,inits,parameters,"robust_cjs.bug",n.chains=nc,n.thin=nt,n.iter=ni,n.burnin=nb)
date() # to tell how long it ran
# took 3 minutes


#summarize posteriors
print(robust.cjs,digits=3)  
# does ok


library(mcmcplots)
library(MCMCvis)
# basic plots to check convergence and credible intervals,
mcmcplot(robust.cjs,
         parms =  c("mean.c", 
                    "mean.p",
                    "mean.phi")
)


# MCMC vis package   
MCMCplot(robust.cjs,
         xlim = c(-1, 1),
         ref_ovl = TRUE,
         params = c("mean.c", 
                    "mean.p",
                    "mean.phi")
         
)

library(MCMCvis)

out$BUGSoutput$sims.list$D

MCMCplot(robust.cjs)
mcmcplots::mcmcplot(robust.cjs)

#CJS SCR Model------------------------------------------------------
#NOT ROBUST

sink("scr_cjs.bug")
cat("					######<--------------------- uncomment 
model{
	
###############Priors and constraints
mean.phi ~ dunif(0, 1) # prior for phi
#mean.p ~ dunif(0, 1)   # prior for p
alpha0 ~ dnorm(0,0.1) 
logit(mean.p) <- alpha0
alpha1 ~ dunif(-50,50)


    
for(i in 1:nind){
	
	for(m in f[i]:n.primary.occasions){  ### for p need every time for phi need -1
		# phi has only 2 dimensions [indiv, and primary occasions]
    phi[i,m] <- mean.phi   # could specify covariates here
    # p and c have 3 dimensions [indiv, secondary, primary]		
    
    for(d in 1:max(n.secondary.occasions)){
      p[i, d, m] <- mean.p  # could specify temporal or individual covariates here
	   } #d for days
    } #m for months
	} #i for individual

#############Likelihood 		

for(i in 1:nind){
	#define latent state at first capture 
  # dimensions [individual, primary session (month)]
	z[i,f[i]] <- 1		# z is true (latent) state alive or dead, know alive at first capture
	
	for(m in  (f[i]+1):n.primary.occasions){  
		
	  #state process				# alive or dead
	  z[i, m] ~ dbern(mu1[i, m]) 		#mu1 is probability alive
		mu1[i, m] <- phi[i, m] * z[i, m-1] # this assures that animals stay dead
	} # m
	
	for(m in prim.occasions[f.prim[i]:length(prim.occasions)] ) { 
	  #observation process			# caught or not
	  #first secondary occasion within a primary occasion:
	  y[i, 1, m] ~ dbern(p.eff[i, 1, m])
	  p.eff[i, 1, m] <- z[i, m] * p[i, 1, m]   
	  
	  #loop over rest of secondary occasions per primary
	  for(d in 2:n.secondary.occasions[m]){ ## add index for number of secondary occasions if not consistent
		  y[i, d, m] ~ dbern(p.eff[i, d, m]) 		# p.eff is prob of capture
		  p.eff[i, d, m] <- z[i, m] * ifelse(sum(y[i, 1:(d-1), m])==0, p[i, d, m], c[i, d, m])	
      # capture prob= (p if not caught previously that session or c if was caught that session) 
      # multiply by if it was alive (so can't capture animal that's not alive)
		 
	   } #d
		} #m
	} #i
}
",fill=TRUE)  #####<----------------uncomment this
sink()



#Run Model----------------------------------------------------------------------
library(jagsUI)

date()
## Call JAGS from R
robust.cjs=jags(bugs.data,inits,parameters,"robust_cjs.bug",n.chains=nc,n.thin=nt,n.iter=ni,n.burnin=nb)
date() # to tell how long it ran
# took 3 minutes


#summarize posteriors
print(robust.cjs,digits=3)  
# does ok


library(mcmcplots)
library(MCMCvis)
# basic plots to check convergence and credible intervals,
mcmcplot(robust.cjs,
         parms =  c("mean.c", 
                    "mean.p",
                    "mean.phi")
)


# MCMC vis package   
MCMCplot(robust.cjs,
         xlim = c(-1, 1),
         ref_ovl = TRUE,
         params = c("mean.c", 
                    "mean.p",
                    "mean.phi")
         
)

library(MCMCvis)

out$BUGSoutput$sims.list$D

MCMCplot(robust.cjs)
mcmcplots::mcmcplot(robust.cjs)
