e2dist <- function (x, y)
{ i <- sort(rep(1:nrow(y), nrow(x)))
dvec <- sqrt((x[, 1] - y[i, 1])^2 + (x[, 2] - y[i, 2])^2)
matrix(dvec, nrow = nrow(x), ncol = nrow(y), byrow = F)
}

#observation array 
M = 150
N = 40
T = 5
K = 5
sigP = 0.5 #within year scale parameter for encounter probability (how much movement)
sigS = 0.5 #alpha 1 base number
phi0 = 0.75
lam0 =0.9

#Built objects and matricies
nreps<- K #number of secondary captures
lam0<-rep(lam0,T) #baseline encounter rate, same for all time periods
phi<-rep(phi0,T) #survival, same for all time periods
pmat<-lam<-list() #p = prob of cap, lam = encounter rate adjusted


#set up the basic trap array in a 10x10 grid
gridx<-seq(1,10,1)
grid<-as.matrix(expand.grid(gridx,gridx))
J=dim(grid)[1]
ntraps = J

#set up space dimensions S
buffer<-2
xl<-min(grid[,1] - buffer)
xu<-max(grid[,1] + buffer)
yl<-min(grid[,2] - buffer)
yu<-max(grid[,2] + buffer)

#Set up individual centroids t=1 
sx<-sy<-sin<-matrix(NA, nrow=M, ncol=T) #rows = ind, col = T
sx[,1]<-runif(M,xl,xu)
sy[,1]<-runif(M,yl,yu)

#set up for initial values for t=1
gamma<-NULL #recruitment probability
gamma[1]<- N/M #recruitment probability at each time step 
z<-r<-matrix(0,nrow=M,ncol=T) #for each individual in time T
r[,1]<-rbinom(M,1,gamma[1])#r = matrix of who get recruited at each time tep
z[,1]<-r[,1]  #z is if its alive at each time step 

#For All other Capture Events after the first one
for (t in 2:T){
  sx[,t]<-runif(M,xl,xu) #generate the x coords
  sy[,t]<-runif(M,yl,yu) #generate the y coords
  surv<-rbinom(M,1,z[,t-1]*phi[t])#look at each individual and see if it survives to the next period
  #generates a list of who survived at each time step

  al<-apply(matrix(z[,1:(t-1)],nrow=M,byrow=FALSE),1,sum)>0 #indicator variable for being alive last time period
  #a1 = 1, true availible to be recruited bc dead last time z[t-1] = 0
  #a1 = 0,false  not availible to be recruited bc alive z[t-1] = 1
  
  #corr acitvity centers for guys alive before
  sx[al,t]<-rtnorm(sum(al),sx[al,t-1],sigS,#movement parameter,
                   lower=xl, upper=xu) #look 
  sy[al,t]<-rtnorm(sum(al),sy[al,t-1],sigS, lower=yl, upper=yu)
 
   idx<- 1- as.numeric(al) #make true/false A1 into binary 1/0 variable
  
   #Recalculate gamma for each Time period 
   gamma[t]<- (N - sum(surv))/sum(idx) #recalculate recruitment 
  
   if (gamma[t]<0) gamma[t]<-0 #what to do if gamma is 0
  r[,t]<- rbinom(M,idx,gamma[t])
  z[,t]<- surv + r[,t] #those who survives plus those recruited 
}

#calculate new S and capture probability per time period 
for (t in 1:T){
  S<-cbind(sx[,t], sy[,t])
  dmat<-e2dist(S,grid)
  psi<- exp(-(1/(2*sigP*sigP))*dmat*dmat) #gaussian random walk 
  lam[[t]]<- lam0[t]*psi
  pmat[[t]]<- 1- exp(-lam[[t]]) #capture probability 
}

#create observation array 
 y<-array(0,dim=c(M,ntraps,T))
    for(t in 1:T){
      yfull<-array(0,dim=c(M,ntraps,K))
      for(i in 1:M){
        for (k in 1:K) {
          yfull[i,1:ntraps,k]<-rbinom(ntraps,1, pmat[[t]][i,]*z[i,t] )
        }
      }
      y[, 1:ntraps,t]<-apply(yfull,1:2,sum)
    }
 #remove individuals not seen 
    ycapt=y[which(rowSums(y[,,])>0),,]

    #calcuate the period of first capture
    first=NULL
    capsums<-apply(ycapt, 3, rowSums)
    for(i in 1:dim(ycapt)[1]){
      a=which(capsums[i,] >0)
      first[i] = min(a)
    }
    
    #function to creat initial S for each year
    Sin<-function(first=first,T=T, M=M, xl=xl, xu=xu, yl=yl, yu=yu,ntot=ntot){
      SX<-SY<-matrix(NA, nrow=ntot, ncol=T)
      for(i in 1:ntot){
        for (t in simdat$first[i]:T){
          SX[i,t]<-runif(1, xl, xu)
          SY[i,t]<-runif(1, yl, yu)
          traps<-which(simdat$y[i,,t]>0)
          if(length(traps)>0){
            SX[i,t]<- mean(grid[traps,1])
            SY[i,t]<- mean(grid[traps,2])
          }
        }
      }
      return(list(SX, SY))
    }
    #set up intial values for z
    zin<-matrix(1, dim(simdat$y)[1], T)
    for(i in 1:dim(simdat$y)[1]){
      for(t in 1:simdat$first[i]){
        zin[i,t]<-NA
      }}
 #statespace is 4x, as is data generation
    data<-list(y=y, first=first, M=dim(y)[1],T=5,xl=xl, xu=xu, yl=yl, yu=yu,K
               =K, J=J, trapmat=as.matrix(grid))
    inits = function() {list(phi=runif(1, .5, 1),sigP=runif(1,1,2), lam0=runif(1), z=zin, sigS=runif(1,0.1,3),
                             SX=Sin(first=simdat$first,T=T, xl=xl, xu=xu, yl=yl, yu=yu, ntot=dim(simdat$y)[1])[[1]],
                             SY=Sin(first=simdat$first,T=T, xl=xl, xu=xu, yl=yl, yu=yu, ntot=dim(simdat$y)[1])[[2]]
    ) }
    params<-c("phi","sigP", "sigS", "lam0")
    cjs.mod2 <- textConnection(cjs.mod)
    mod<-jags.model(cjs.mod, data, inits, n.chains=3, n.adapt=100)
    mod.out <- jags.model("cjs.txt", data, inits, n.chains=3)
    out<-coda.samples(mod.out, params, n.iter=100, thin=3)
    out
    close(cjs.mod)
    library(MCMCvis)
    library(mcmcplots)    