library(AHMbook)
library(jagsUI) 
library(MCMCvis)
library(R2jags)

#Functions-------------------------------------------------------------------------

reality <- function (N, #Number of indididuals
                     trap_no, #Number of traps
                     buffer, #buffer around S 
                     alpha0, #HR alpha0
                     sigma #sigma 
                     ) {
 
#Keep it all the same
  set.seed(1234)
  #create traps and locations
  traplocs<- cbind(sort(rep(1:trap_no,trap_no)),rep(1:trap_no,trap_no))
  Dmat<-e2dist(traplocs,traplocs)
  ntraps<-nrow(traplocs)
  
  #Define S as larger than trapping grid 
  buffer<-2
  Xl<-min(traplocs[,1] - buffer)
  Xu<-max(traplocs[,1] + buffer)
  Yl<-min(traplocs[,2] - buffer)
  Yu<-max(traplocs[,2] + buffer)
  
  #Generate Centroids
  sx<-runif(N,Xl,Xu) #centroid x 
  sy<-runif(N,Yl,Yu) #centroid y
  S<-cbind(sx,sy) #Simulated S with centroids
  
  #Set HR equation 
  alpha0<- alpha0
  sigma<- sigma
  alpha1<- 1/(2*sigma*sigma)
  
  #Generate trapping probability given centroid location per individual
  D<- e2dist(S,traplocs)# how far is each individual from each trap?
  
  probcap<-plogis(alpha0)*exp(-alpha1*D*D)
  
  
  plot(S, main = "S", ylim = c(0,Yu), xlim = c(0, Xu))
  points(traplocs, col = "Red")
  
  
 return(list(traplocs=traplocs,xlim=c(Xl,Xu),ylim=c(Yl,Yu),N=N,alpha0=alpha0,alpha1=alpha1,sigma=sigma, S=S, ntrap = ntraps, probcap=probcap, D=D))

}
observation <-
  function(output, #list from previous function
           cap.no #number of trapping sessions 
  ){
    set.seed(1234)
    cap.no = cap.no
    ntraps = output$ntrap
    probcap = output$probcap
    N = output$N
    
    #Create Capture History
    #SCR Version
    #  Y<-matrix(NA,nrow=N,ncol=ntraps)
    # for(i in 1:nrow(Y)){
    #  Y[i,]<-rbinom(ntraps,cap.no,probcap[i,])
    
    #My version
    Y <- matrix(NA,nrow=N,ncol=ntraps)
    for(i in 1:nrow(Y)){ #individual
      for(j in 1:ncol(Y)){ #trap
        Y[i,j] <- rbinom(1, cap.no, probcap[i,j]) 
        
      }
    }
    
    
    #Remove rows of data with all 0s 
    Y =  Y[rowSums(abs(Y)) != 0,]
    #Spit out matrix
    Y
  }



output <- reality(N = 30, trap_no =10, alpha0 = -0.3, sigma = 3)

observed  <- observation(output,10)
observed

#Plot Si versus Trap Location--------------------------------------------------------
library(ggplot2)

color <- function(df){ #output collection of objects produced initially 
  grid <- data.frame(df$traplocs)
  S <- df$S
  col <- df$probcap
  
  for(i in 1:nrow(col)){
  p <- ggplot(grid, aes(x = X1, y = X2, color= col[i,])) +
    geom_point(size = 3) +
    scale_color_continuous()+
    theme_minimal()+
    annotate("point", size = 3, x = S[i,1], y = S[i,2], colour = "red")
  plot(p)
  }
 
}

color(output)

#Graph Outputs-----------------------------------------------------------------

#Gaussian Output Function
#
plot(output$D,output$probcap) # y= p, x = distance from trap 


#Known N Model-----------------------------------------------------------------
## Grab encounter histories from simulated data list 
y <- observed 

## Grab the trap locations 
traplocs <- output$traplocs

nind <- nrow(y) 
X <- output$traplocs 
J <- nrow(X) #number of traps 
 # number of traps 

xlim <- output$xlim 
ylim <- output$ylim

### Starting values for activity centers,s 

sst <- cbind(runif(nind,xlim[1],xlim[2]),runif(nind,ylim[1],ylim[2])) 

for(i in 1:nind){ #check this here 
  if(sum(y[i,])==0) next
  sst[i,1] <- mean( X[y[i,]>0,1] ) 
  sst[i,2] <- mean( X[y[i,]>0,2] ) }


sink("KnownN.bug") 
cat(" 
     model{ 
     alpha0 ~ dnorm(0,0.1) 
     logit(p0) <- alpha0
     alpha1 ~ dunif(-50,50)
   
     
     
     for(i in 1:N){
     # note N here -- N is KNOWN in this example 
     s[i,1] ~ dunif(xlim[1],xlim[2]) 
     s[i,2] ~ dunif(ylim[1],ylim[2]) 
     
     for(j in 1:J){ 
     d[i,j] <- pow(pow(s[i,1]-X[j,1],2) + pow(s[i,2]-X[j,2],2),0.5)
     y[i,j] ~ dbin(p[i,j],K)
     p[i,j] <- p0*exp(alpha1*d[i,j]*d[i,j]) 
     }
     } #
    
    #sigma <- sqrt(1/(2*alpha1))
     } #model
     ", fill = TRUE)
     sink()

     
data <- list (y=y, X=X, K=K, N=nind, J=J, xlim=xlim, ylim=ylim) 
inits <-  function() {list (alpha0=rnorm(1,-4,.4), alpha1=runif(1,1,2), s=sst, z =z) } #Adapt to Jags 
inits
     
library(jagsUI) 
library(R2jags)
     
parameters <- c("alpha0","alpha1","sigma") 
     
out <- jags(data = data, inits = NULL, parameters = parameters, "KnownN.bug", n.thin=1, n.chains=3, n.burnin=1000,n.iter=10000)
out
   
mcmcplots::mcmcplot(out)  
#Unknown N; Augmentation------------------------------------------------------ 

y <-observed

y

     
space <- output$ylim*output$xlim
nind <- nrow(y) 
X <- output$traplocs 
J <- nrow(X) #number of traps 
 # number of traps 
#K <- 20
xlim <- output$xlim 
ylim <- output$ylim

M = 50 #Augmented population
y = rbind(y, matrix(0, nrow=M-nind, ncol = ncol(y)))    #Augmented capture history

z <-c(rep(1,nind),rep(NA,M-nind))
z
### Starting values for activity centers,s 
     
sst <- cbind(runif(M,xlim[1],xlim[2]),runif(M,ylim[1],ylim[2])) 
for(i in 1:nind){ #check this here 
       if(sum(y[i,])==0) next
       sst[i,1] <- mean( X[y[i,]>0,1] ) 
       sst[i,2] <- mean( X[y[i,]>0,2] ) }

          
     
sink("DataAugmentation.bug") 
     cat(" 
     model{ 
     alpha0 ~ dunif(-50,50) 
     logit(p0) <- alpha0 
     alpha1 ~ dunif(-50,50)
     #sigma <- sqrt(1/(2^alpha1))
     psi ~ dunif(0,1)
     
     for(i in 1:M) {
     z[i] ~ dbern(psi)
     s[i,1] ~ dunif(xlim[1],xlim[2]) 
     s[i,2] ~ dunif(ylim[1],ylim[2]) 
     
     for(j in 1:J){ 
     d[i,j] <- pow(pow(s[i,1]-X[j,1],2) + pow(s[i,2]-X[j,2],2),0.5)
     y[i,j] ~ dbin((p.eff[i,j]),K)
     p.eff[i,j] <- z[i]*p[i,j]
     p[i,j] <- p0*exp(alpha1*d[i,j]*d[i,j]) 


     }
     } #observations
     N <- sum(z[1:M])
     D <- N/144  #(xlim[2]*ylim[2])
     } #model
     ", fill = TRUE)
     sink()


 data <- list (y=y, X=X, K=K, J=J, xlim=xlim, ylim=ylim, M=M, z =z) 
 data
 inits <- function(){
   list(z=z, s=sst, alpha0=rnorm(1,-4,.4), alpha1=runif(1,1,2), psi= runif(1))
 }
inits
 parameters <- c("alpha0","alpha1", "psi", "N", "D") 
 
 out <- jags(data = data, inits = NULL, parameters = parameters, "DataAugmentation.bug", n.thin=1, n.chains=3, n.burnin=1000,n.iter=10000)
 out
 
 library(MCMCvis)
 
 out$BUGSoutput$sims.list$D

 MCMCplot(out)
 mcmcplots::mcmcplot(out)

#Value Testing------------------------------------------------------------------

k_testing <- function(N, #size of population
                      number) #highest value for K
                      { 
#Select vector with values and which parameter
N <- N
k <- seq(1,number, by = 1)
#Build final value matrix
values <- data.frame(matrix(NA, nrow = length(k), ncol = 4))
colnames(values) <- c("k", "N", "Nob", "percent")

values$k <- k
values$N <- N 

output <- reality(N = N, trap_no =10, alpha0 = -0.5, sigma = 3)
#simulate on that, keeping everything the same 
for (i in 1:length(k)){
  #Simulate Data
  #Currently all the same reality, different # of trapping occasions
    observed  <- observation(output,k[i])
  #Run Model
    y <-observed #capture history
    nind <- nrow(y) #number of individuals
    X <- output$traplocs #trap locations
    J <- nrow(X) #number of traps 
    xlim <- output$xlim #Boundary of S
    ylim <- output$ylim #Boundary of S
    M = 50 #Augmented population
    y = rbind(y, matrix(0, nrow=M-nind, ncol = ncol(y)))    #Augmented capture history
    z <- c(rep(1,nind),rep(NA,M-nind))
    ### Starting values for activity centers,s 
    sst <- cbind(runif(M,xlim[1],xlim[2]),runif(M,ylim[1],ylim[2])) 
    for(j in 1:nind){ #check this here 
        if(sum(y[j,])==0) next
        sst[j,1] <- mean( X[y[j,]>0,1] ) 
        sst[j,2] <- mean( X[y[j,]>0,2] ) } 
   #Load Data for Jags
     data <- list (y=y, X=X, K=k[i], J=J, xlim=xlim, ylim=ylim, M=M, z=z) 
   #Identify Parameters
    parameters <- c("alpha0","alpha1", "psi", "N", "D") 
   #Run Model 
    out <- jags(data = data, inits = NULL, parameters = parameters, "DataAugmentation.bug", n.thin=1, n.chains=3, n.burnin=1000,n.iter=2000)
  #Put estimated N in values matrix
    values[i,3] = out$BUGSoutput$mean$N 
print(k[i])
}

#Calculate Missed N
values$percent <- values$N/values$Nob
#Plot Matrix 
p <- plot(values$k,values$percent)
#output
p
}

#Actually Run it 
Sys.time()
output_test <- k_testing(30,100)

Sys.time()

alpha0_test <- function(N, #size of population
                        low, #lowest value of alpha0
                        high, #highest value of alpha0
                        K) #number of captures
{ 
  #Select vector with values and which parameter
  N <- N
  k <- K
  a0 <- seq(low,high, by = 0.25)
  #Build final value matrix
  values <- data.frame(matrix(NA, nrow = length(a0), ncol = 4))
  colnames(values) <- c("a0", "N", "Nob", "percent")
  
  values$a0 <- a0
  values$N <- N 
  

  #simulate on that, keeping everything the same 
  for (i in 1:length(a0)){
    #Simulate Data
    #Currently all the same reality, different # of trapping occasions
    output <- reality(N = N, trap_no =10, alpha0 = a0[i], sigma = 3)
    observed  <- observation(output,k)
    #Run Model
    y <-observed #capture history
    nind <- nrow(y) #number of individuals
    X <- output$traplocs #trap locations
    J <- nrow(X) #number of traps 
    xlim <- output$xlim #Boundary of S
    ylim <- output$ylim #Boundary of S
    M = 50 #Augmented population
    y = rbind(y, matrix(0, nrow=M-nind, ncol = ncol(y)))    #Augmented capture history
    z <- c(rep(1,nind),rep(NA,M-nind))
    ### Starting values for activity centers,s 
    sst <- cbind(runif(M,xlim[1],xlim[2]),runif(M,ylim[1],ylim[2])) 
    for(j in 1:nind){ #check this here 
      if(sum(y[j,])==0) next
      sst[j,1] <- mean( X[y[j,]>0,1] ) 
      sst[j,2] <- mean( X[y[j,]>0,2] ) } 
    #Load Data for Jags
    data <- list (y=y, X=X, K=k, J=J, xlim=xlim, ylim=ylim, M=M, z=z) 
    #Identify Parameters
    parameters <- c("alpha0","alpha1","sigma", "psi", "N", "D") 
    #Run Model 
    out <- jags(data = data, inits = NULL, parameters = parameters, "DataAugmentation.bug", n.thin=1, n.chains=3, n.burnin=1000,n.iter=2000)
    #Put estimated N in values matrix
    values[i,3] = out$BUGSoutput$mean$N 
    print(a0[i])
  }
  
  #Calculate Missed N
  values$percent <- values$N/values$Nob
  #Plot Matrix 
  q <- plot(values$a0,values$percent)
  #plot detection function
  p <- plot(output$D,output$probcap) # y= p, x = distance from trap
  
  p
  q
}
Sys.time()
alpha0_test(30, -1, 5, 5)
Sys.time()

sigma_test <- function(N, #size of population
                        low, #lowest value of sigma (must be 2 or larger)
                        high, #highest value of sigma
                        K) #number of captures 
{ 
  #K = 3
  #low = 3
  #high = 10
 # N = 30
  #Select vector with values and which parameter
  N <- N
  k <- K
  sigma <- seq(low,high, by = 2)
  #Build final value matrix
  values <- data.frame(matrix(NA, nrow = length(sigma), ncol = 4))
  colnames(values) <- c("sigma", "N", "Nob", "percent")
  
  values$sigma <- sigma
  values$N <- N 
  
  
  #simulate on that, keeping everything the same 
  for (i in 1:length(sigma)){
    #Simulate Data

    #Currently all the same reality, different # of trapping occasions
    output <- reality(N = N, trap_no=10, alpha0 = -0.5, sigma = sigma[i])
    observed  <- observation(output,k)
    #Run Model
    y <-observed #capture history
    nind <- nrow(y) #number of individuals
    X <- output$traplocs #trap locations
    J <- nrow(X) #number of traps 
    xlim <- output$xlim #Boundary of S
    ylim <- output$ylim #Boundary of S
    M = 50 #Augmented population
    y = rbind(y, matrix(0, nrow=M-nind, ncol = ncol(y)))    #Augmented capture history
    z <- c(rep(1,nind),rep(NA,M-nind))
    ### Starting values for activity centers,s 
    sst <- cbind(runif(M,xlim[1],xlim[2]),runif(M,ylim[1],ylim[2])) 
    for(j in 1:nind){ #check this here 
      if(sum(y[j,])==0) next
      sst[j,1] <- mean( X[y[j,]>0,1] ) 
      sst[j,2] <- mean( X[y[j,]>0,2] ) } 
    #Load Data for Jags
    data <- list (y=y, X=X, K=k, J=J, xlim=xlim, ylim=ylim, M=M, z=z) 
    #Identify Parameters
    parameters <- c("alpha0","alpha1","sigma", "psi", "N", "D") 
    #Run Model 
    out <- jags(data = data, inits = NULL, parameters = parameters, "DataAugmentation.bug", n.thin=1, n.chains=3, n.burnin=1000,n.iter=2000)
    #Put estimated N in values matrix
    values[i,3] = out$BUGSoutput$mean$N 
    print(sigma[i])
  }
  
  #Calculate Missed N
  values$percent <- values$N/values$Nob
  #Plot Matrix 
  q <- plot(values$sigma,values$percent)
  #plot detection function
  p <- plot(output$D,output$probcap) # y= p, x = distance from trap
  
  p
  
}

Sys.time()
sigma_testing <- sigma_test(30, 2, 50, 5)
Sys.time()

