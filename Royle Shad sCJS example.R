#install required packages
library(rjags)
library(reshape)
data(Ch16shaddata)

##data setup and organization
# Input and format data matrix for spatial CJS
shad <- Ch16shaddata$shad
melted.rkm <- melt(shad, id=c("TagID","RKM")) 
tagid.week <- cast(melted.rkm, TagID ~ RKM ~ value, fill=0, length)
tagid.week[,,13]

#data for non-spatial CJS
hold=as.matrix(table(shad$TagID, shad$Week))
hold[hold > 1] <- 1
y <- hold

first<-Ch16shaddata$first


# Constants:
M <- 315       # Number of individuals
T <- 12     # Number of periods (weeks)

#set up trap locations and state space
nantenna <- 7  # weir, 6 antennas
antenna.loc <- c(3.7, 7.7, 13.4, 45.3, 56.4, 72.0, 77.0)  # antenna locations
xl <- 0     # lower boundary, river mouth
xu <- 82    # upper boundary, Atkinson Mill Dam

######################################################################################################
###### Non-Spatial Cormack Jolly Seber   NS-CJS#######################################################


  sink("ModelNSCJS.txt")
  cat("

model { #model
# Priors

phi ~ dunif(0,1)   # Survival (constant)

for(t in 1:T){
p[t] ~ dunif(0, 1)    # detection (constant)
}
    
for (i in 1:M){ #m    

   z[i,first[i]] ~ dbern(1)  # Known to be alive at entry into study area

  for (t in (first[i]+1):T) { #t
       #State process--is it alive? 
           z[i,t] ~ dbern(phiUP[i,t]) #PhiUP = Mu1, probability of being alive
       phiUP[i,t] <- z[i,t-1]*phi #Ensures that the animal stays dead 
       #Observation process--was it captured?
           y[i,t] ~ dbern(tmp[i,t])
         tmp[i,t] <- z[i,t]*p[t] #mu2 = tmp,  Angie writes as p.eff

     

         
   } # t
  } # m
} #model


", fill = TRUE)
  sink()
  
  
  ######################################################
  
  #Set up data input
  data<-list(y=y, first=first, M=M, T=T)
  
  z=matrix(NA, M, T)
  for(i in 1:M){ 
    for(t in first[i]:T){ 
      z[i,t] <-1
    } 
  }
  
  #Set initial values
  inits =  function() {list(z=z, phi=runif(1,0,1), p=runif(12,0,1)) }
  
  # Parameters to follow
  parameters <- c("phi", "p")
  
  mod.out <- jags.model("ModelNSCJS.txt", data, inits, n.chains=n.chains)
  out <- coda.samples(mod.out,  parameters, n.iter=10000)                                                  
  
} #end if

################################################################################
#Spatial CJS
################################################################################
sink("SpatialCJS")
cat("

model {
# Priors
sigma ~ dunif(0,80)
sigma2 <- sigma*sigma
phi ~ dunif(0, 1)   # Survival (constant across time)
tauv~dunif(0, 30)
tau<-dunif(0,30)

for(t in 1:T){
p0[t] ~ dgamma(0.1, 0.1)
}

for (i in 1:M){
  z[i,first[i]] <- 1
  S[i,first[i]] ~ dunif(0,50)

for(j in 1:nantenna) {
  D2[i,j,first[i]] <- pow(S[i,first[i]]-antenna.loc[j], 2)
       p[i,j,first[i]]<-  p0[first[i]]*exp(- D2[i,j,first[i]]/(2*sigma2))
       tmp[i,j,first[i]] <- p[i,j,first[i]]
         y[i,j,first[i]] ~ dbern(tmp[i,j,first[i]])
      }

   for (t in first[i]+1:T) {
          S[i,t] ~ dunif(xl, xu)
         for(j in 1:nantenna) {
                D2[i,j,t] <- pow(S[i,t]-antenna.loc[j], 2)
               p[i,j,t] <- p0[t] * exp(-D2[i,j,t]/(2*sigma2))
               tmp[i,j,t] <- z[i,t]*p[i,j,t]
                 y[i,j,t] ~ dbern(tmp[i,j,t])
 }
    phiUP[i,t] <-  z[i,t-1]*phi
       z[i,t] ~ dbern(phiUP[i,t])
}
}
}
", fill = TRUE)
sink()

######################################################

#Set up data input

data<-list(y=tagid.week, first=first, M=M, T=T, xl=xl, xu=xu, nantenna=nantenna, antenna.loc=antenna.loc)

z=matrix(NA, M, T)
St=matrix(NA, M, T)

for(i in 1:M){ 
  for(t in (first[i]+1):T){ 
    z[i,t] <-1
    St[i,t] <- (sum(tagid.week[i,,t] * antenna.loc)/sum(tagid.week[i,,t])) + 3
    if(St[i,t] == "NaN") St[i,t] <- St[i,t-1]
  } 
}
z
St

#Set initial values
inits =  function() {list(z=z, phi=runif(1,0,1), p0=runif(12,1,3), tauv=runif(1,0,20)) }

# Parameters to follow
parameters <- c("sigma", "phi", "p0")

mod.out <- jags.model("SpatialCJS", data, inits, n.chains=3)

 n.chains = 3
 