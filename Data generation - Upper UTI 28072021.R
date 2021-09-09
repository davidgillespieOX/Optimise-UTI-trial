#For optimise UTI trial 08042021

# Load relevant libraries
library(boot)

# Set seed
set.seed(1)

# Load Simulation Scenarios
source("C:/Users/wmddg1/Documents/University of Oxford/Grant applications/HTA antibiotic duration/Stage 2/Simulations for sample size/Upper UTI/Simulation Scenarios UTI - Upper UTI 28072021.R")

# Set working directory
setwd("C:/Users/wmddg1/Documents/University of Oxford/Grant applications/HTA antibiotic duration/Stage 2/Simulations for sample size/Upper UTI")

# Load Simulation Parameters
load("C:/Users/wmddg1/Documents/University of Oxford/Grant applications/HTA antibiotic duration/Stage 2/Simulations for sample size/Upper UTI/Simulation Parameters Upper UTI 28072021.RData")
durations<-c(4, 6, 8, 10, 12, 14)
n.arms<-6
n.obj<-n<-M.boot<-500


for (s in 1:n.scen) {
  
  data<-array(NA,c(n,2,n.sim))
  formula.1<-get(paste("formula",scenarios[s],".4", sep=""))

  
  #Calculate expected outcome given scenario
  y.dur<-formula.1(x.dur)
  max.p<-formula.1(max.dur)
  
  for (i in 1:n.sim) {
    
    
    # Durations in long format 
    durlong<-sample(durations,n,rep=TRUE)
    
    # Calculate probabilities
    pro<-formula.1(durlong)
    
    # Generate events from binomial distribution
    y<-rbinom(n,1,pro)
    
    data[,1,i]<-durlong
    data[,2,i]<-y
    
  }
  assign(paste("Scenario",s,sep=""),data)
}

whole.data<-mget(paste("Scenario", 1:n.scen, sep=""))
save(whole.data, file="Simulated Data Upper UTI 28072021.RData")