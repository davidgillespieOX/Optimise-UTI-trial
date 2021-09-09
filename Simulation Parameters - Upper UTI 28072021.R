# Clear R environment
rm(list=ls())

# Set Working directory
setwd("C:/Users/wmddg1/Documents/University of Oxford/Grant applications/HTA antibiotic duration/Stage 2/Simulations for sample size/Upper UTI")

# Parameter values

n.sim<-20
n.obj<-n<-500
n.arms<-6
#npergroup<-ceiling(n.obj/n.arms)
#n<-npergroup*n.arms
scenarios<-c(1:8)
n.scen<-length(scenarios)
min.dur<-4
max.dur<-14
durations<-seq(min.dur,max.dur,length.out = n.arms) # We choose equidistant duration arms, between the minimum and the maximum
poss.durations<-seq(min.dur,max.dur-1)
all.durations<-seq(min.dur,max.dur)
n.dur<-length(all.durations)
#durlong<-rep(durations, each=npergroup)
x.dur<-seq(min.dur,max.dur, length.out = 100)
alpha<-0.025
M.boot<-500

# Save Parameters

save.image(file="Simulation Parameters Upper UTI 28072021.RData")
