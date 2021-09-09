# Load libraries
library(splines)
library(boot)
library(flux)
library(mfp)

#Replication of Table 2 from 2018 Quartagno paper

# Set seed
set.seed(1)

#Load Simulation Scenarios
source("C:/Users/wmddg1/Documents/University of Oxford/Grant applications/HTA antibiotic duration/Stage 2/Simulations for sample size/Upper UTI/Simulation Scenarios UTI - Upper UTI 28072021.R")

# Set working directory
setwd("C:/Users/wmddg1/Documents/University of Oxford/Grant applications/HTA antibiotic duration/Stage 2/Simulations for sample size/Upper UTI")

#Some parameters initializations
n.sim<-1000
n.obj<-500
n.arms<-6
scenarios<-c(1:8)
n.scen<-length(scenarios)
min.dur<-4
max.dur<-14
x.dur<-seq(min.dur,max.dur, length.out = 100)

# Some matrices initializations:
max.err.fp<-matrix(NA,n.sim,n.scen)
where.max.err.fp<-matrix(NA,n.sim,n.scen)
diff.auc.fp<-matrix(NA,n.sim,n.scen)
reldiff.auc.fp<-matrix(NA,n.sim,n.scen)
cover.auc.fp<-matrix(NA,n.sim,n.scen)
list.of.fits<-vector(mode="list", length=n.scen*n.sim)
all.data<-array(NA,dim=c(ceiling(n.obj/n.arms)*n.arms,n.sim,n.scen))

# Start simulations looping over different scenarios
for (s in 1:n.scen) {
  formula.1<-get(paste("formula",scenarios[s],".4", sep=""))
  
  #Calculate y of formula to calculate AUC
  y.dur<-formula.1(x.dur)
  
  #Calculate true area under the curve
  true.auc<-auc(x.dur,y.dur)
  
  # Define formula to calculate root in prediction curve for fractional polynomial model:
  formula.3<-function(x) { 
    return(predict(fit.fp.i,data.frame(durlong=x),type="resp"))
  }
  
  
  # Run n.sim simulations:
  for (i in 1:n.sim) {
    
    #Do same analysis for different number of duration arms
    
    
    # We choose equidistant duration arms, between the minimum and the maximum
    durations<-seq(min.dur,max.dur,length.out = n.arms)
    
    # number of patients in each duration group, aiming for a total of n.obj patients
    npergroup<-ceiling(n.obj/n.arms)
    
    # total number of patients
    n<-npergroup*n.arms
    
    # Durations in long format 
    durlong<-rep(durations, each=npergroup)
    
    # Calculate probabilities
    pro<-formula.1(durlong)
    # Generate events from binomial distribution
    y<-rbinom(n,1,pro)
    
    # Fit Fractional Polynomial regression model
    fit.fp.i<-mfp(y~fp(durlong), family="binomial")
    list.of.fits[[(s-1)*1000+i]]<-fit.fp.i
    
    #calculate area between true and estimated curve:
    
    y.dur.est<-formula.3(x.dur)
    y.dur.ses<-predict(fit.fp.i,data.frame(durlong=x.dur),type="resp", se.fit=T)$se.fit
    y.fit<-formula.3(x.dur)-y.dur
    est.auc<-auc(x.dur,abs(y.fit))
    max.err<-max(abs(y.fit))
    where.max.err<-which.max(abs(y.fit))
    cov<-0
    for (j in 1:100) {
      if (y.dur.est[j]-qnorm(0.975)*y.dur.ses[j]<y.dur[j]&y.dur.est[j]+qnorm(0.975)*y.dur.ses[j]>y.dur[j]) {
        cov<-cov+1
      }
    }
    
    #results
    max.err.fp[i,s]<-max.err
    where.max.err.fp[i,s]<-x.dur[where.max.err]
    diff.auc.fp[i,s]<-est.auc
    reldiff.auc.fp[i,s]<-est.auc/(max.dur-min.dur)
    cover.auc.fp[i,s]<-cov
    all.data[,i,s]<-y
    
    if (i%%10==0) cat("Simulation number ", i, " completed\n")
  }
  cat("Scenario number ", s, " completed\n")
  
}
Table<-matrix(NA,n.scen+1,8)
for (s in 1:n.scen) {
  Table[s,1:5]<-quantile(reldiff.auc.fp[,s], c(0, 0.05, .50, .95, 1))
  Table[s,6:7]<-quantile(max.err.fp[,s], c(.50, .95))
  Table[s,8]<-mean(cover.auc.fp[,s])
}
Table[n.scen+1,1:5]<-quantile(as.vector(reldiff.auc.fp), c(0, 0.05, .50, .95, 1))
Table[n.scen+1,6:7]<-quantile(as.vector(max.err.fp), c(.50, .95))
Table[n.scen+1,8]<-mean(as.vector(cover.auc.fp))

rownames(Table)<-c("Scenario 1", "Scenario 2", "Scenario 3", "Scenario 4", "Scenario 5", "Scenario 6", "Scenario 7", "Scenario 8", "Scenario 9", "Scenario 10", "Scenario 11", "Scenario 12", "Scenario 13", "Scenario 14", "Scenario 15", "Scenario 16", "Overall")
colnames(Table)<-c("Min", "5th Perc", "Median", "95th Perc", "Max", "Median Max Diff", "95th Perc max Diff", "Coverage")
windows() 
par(mfrow=c(4,4))

for (s in 1:n.scen) {
  hist(where.max.err.fp[,s], breaks = c(seq(1,5,0.5)), xlim=c(1,5), main = paste("Scenario ",s))
}
savePlot(filename = "histogram_duration_maxAE", type="png")
dev.off() 

windows()
par(mfrow=c(4,4))

for (s in 1:n.scen) {
  formula.4<-get(paste("formula",scenarios[s],".4", sep=""))
  curve(formula.4(x), xlim=c(min.dur,max.dur), ylim=c(0,1), xlab = "Duration", ylab="% cure", lwd=3, main=paste("Scenario",s))
  worst100<-quantile(reldiff.auc.fp[,s],0.2)
  worst.indexes<-which(reldiff.auc.fp[,s]>worst100)
  for (i in worst.indexes)  curve(predict(list.of.fits[[(s-1)*1000+i]],data.frame(durlong=x),type="resp"),add=T, col="red")
  curve(formula.4(x), xlim=c(min.dur,max.dur),  lwd=3,  add=T)
  
}
savePlot(filename = "worst100Simulations", type="png")

dev.off()
save(Table, all.data, diff.auc.fp, reldiff.auc.fp, where.max.err.fp, max.err.fp,cover.auc.fp, file = "Results.fp.Table1.RData")
save(list.of.fits, file = "listoffits_basecase.RData")

Table
