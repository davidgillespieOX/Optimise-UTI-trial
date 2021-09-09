# An example of duration-response curve:

# Load relevant libraries
library(boot)
library(gamlss)
fp<-gamlss::fp

# Set seed
set.seed(1)

#Load Simulation Scenarios
source("C:/Users/wmddg1/Documents/University of Oxford/Grant applications/HTA antibiotic duration/Stage 2/Simulations for sample size/Upper UTI/Simulation Scenarios UTI - Upper UTI 28072021.R")

# Set working directory
setwd("C:/Users/wmddg1/Documents/University of Oxford/Grant applications/HTA antibiotic duration/Stage 2/Simulations for sample size/Upper UTI")

# Load Simulation Parameters
load("C:/Users/wmddg1/Documents/University of Oxford/Grant applications/HTA antibiotic duration/Stage 2/Simulations for sample size/Upper UTI/Simulation Parameters Upper UTI 28072021.RData")

# Load Simulated Data
load("C:/Users/wmddg1/Documents/University of Oxford/Grant applications/HTA antibiotic duration/Stage 2/Simulations for sample size/Upper UTI/Simulated Data Upper UTI 28072021.RData")

# Set working directory
setwd("C:/Users/wmddg1/Documents/University of Oxford/Grant applications/HTA antibiotic duration/Stage 2/Simulations for sample size/Upper UTI")

# Select scenario 6 and first simulation:
s=4
i=1
y<-whole.data[[s]][,2,i]
durlong<-whole.data[[s]][,1,i]
data.mfp<-data.frame(y, durlong)
fit.i<-gamlss(y~fp(durlong), data=data.mfp, trace=F, family = BI)

# Define formula to predict curve:
predicted<-function(x, fit, datac) { 
  return(predict(object=fit,newdata=data.frame(durlong=x), data=datac, type="resp"))
}
invisible(capture.output(y.dur.est<-predicted(max.dur, fit.i, data.mfp)))

# Run 500 bootstrap replications:

find.min.dur<- function (data.mfp, indices) {
  # Select bootstrap sampel:
  da <- data.mfp[indices,]
  
  # Fit Fractional Polynomial regression model
  fit.i<-gamlss(y~fp(durlong), data=da, trace=F, family = BI)
  
  # Predict duration-response curve and associated pointwise CI
  invisible(capture.output(y.dur.est2<-predicted(x.dur, fit.i, da)))
  invisible(capture.output(y.durations<-predicted(max.dur, fit.i, da)-predicted(poss.durations[1:(n.dur-1)], fit.i, da)))
  
  #Define acceptability curve. Constant at 10% less than control.
  acceptability2<-function(x) {
    return(y.dur.est2[100]-0.1)
  }
  y.accept<-acceptability2(x.dur)
  
  # What is point where predicted lower CI first crosses acceptability curve?
  flag=t=1
  min.duration<-max.dur
  while ((t<length(x.dur))&(flag==1)) {
    if ((y.dur.est2[t]-y.accept)>0) {
      flag=0
      min.duration<-x.dur[t]
    }
    t=t+1
  }
  output<-c(min.duration, y.durations)
  return(output)
}
results<-boot(data.mfp,find.min.dur,M.boot)
if (length(unique(results$t[,1])) > 1) {
  res.ci<-boot.ci(results, conf=1-alpha*2, type="perc", index=1)
  min.duration<- res.ci$perc[4:5]
} else {
  min.duration<-results$t0
}
up.bounds.CI<-NULL
low.bounds.CI<-NULL
for (indx in 2:length(results$t0)) {
  res.ci2<-boot.ci(results, conf=1-alpha*2, type="bca", index=indx)
  up.bounds.CI<-c(up.bounds.CI, res.ci2$bca[5])
  low.bounds.CI<-c(low.bounds.CI, res.ci2$bca[4])
  
}

est.opt.dur<-results$t0[1]
est.opt.cure<-predicted(est.opt.dur, fit.i, data.mfp)

pdf(file="Analysis example for Upper UTI.pdf",width=12,height=7)

par(mfrow=c(1,2))

curve(predicted(x, fit.i, data.mfp), xlim=c(min.dur,max.dur), ylim=c(0,1), xlab = "Duration", ylab="% sustained cure", lwd=3, main="Duration-Response Curve", xaxt="n", yaxt="n")
axis(side=1, at=all.durations, labels=all.durations)
axis(side=1, at=9, labels=9, col.axis="red", col.ticks = "red")
axis(side=2, at=seq(0,1,by=0.2), labels=paste(seq(0,100,by=20),"%"))
abline(h=(y.dur.est-0.1), col="red")
segments(min.duration[1],0.1, min.duration[2], 0.1, lwd=2)
segments(min.duration[2],0.1, 9, 0.1, lwd=1, col="red", lty=2)
points(9,0.1,col="red", pch=8)
segments(9,0.1, 9, -0.05, lwd=1, col="red", lty=2)
segments(min.duration[1],0.1-0.01,min.duration[1],0.1+0.01, lwd=2)
segments(min.duration[2],0.1-0.01,min.duration[2],0.1+0.01, lwd=2)
segments(est.opt.dur,0.1-0.005,est.opt.dur,0.1+0.005, lwd=2)
segments(est.opt.dur,est.opt.cure,est.opt.dur,0.1+0.005, lwd=1, col="grey", lty=2)
text(11,0.87,labels="Duration-Response Curve", font=2, cex=0.75)
text(12,0.75,labels="Acceptability Frontier", col="red", font=2, cex=0.75)
text(6,0.07,labels="Bootstrap CI for Optimal Duration", font=2, cex=0.75)

curve(rep(0.1,length(x)), xlim=c(min.dur,poss.durations[n.dur-1]), ylim=c(-0.2,0.8),xaxt="n", yaxt="n", xlab = "Duration", ylab="Diff. % sustained cure", lty=1, col="red", main="Difference from longest duration")
axis(side=1, at=all.durations, labels=all.durations)
axis(side=1, at=9, labels=9, col.axis="red", col.ticks = "red")
axis(side=2, at=seq(-0.2,0.6,by=0.2), labels=paste(seq(-20,60,by=20),"%"))
for (d in (1:(n.dur-1))) {
  if (d==6) color.plot<-"red" else color.plot<-"black"
  segments(poss.durations[d]-0.1,low.bounds.CI[d], poss.durations[d]+0.1, low.bounds.CI[d], lwd=2, col=color.plot)
  segments(poss.durations[d]-0.05,results$t0[d+1],poss.durations[d]+0.05,results$t0[d+1], lwd=2, col=color.plot)
  segments(poss.durations[d]-0.1,up.bounds.CI[d],poss.durations[d]+0.1,up.bounds.CI[d], lwd=2, col=color.plot)
  segments(poss.durations[d],up.bounds.CI[d],poss.durations[d],low.bounds.CI[d], lwd=2, col=color.plot)
}
segments(poss.durations[6],-0.3,poss.durations[6],low.bounds.CI[6], lwd=1, col="red", lty=2)
text(12,0.14,labels="Acceptability Frontier", col="red", font=2, cex=0.75)
text(7.5,0.38,labels="Bootstrap CI for Difference\n in Efficacy from d=14", font=2, cex=0.75)
text(5,-0.03,labels="Equality Line", font=2, cex=0.75)

abline(h=0, lty=2)
dev.off()
