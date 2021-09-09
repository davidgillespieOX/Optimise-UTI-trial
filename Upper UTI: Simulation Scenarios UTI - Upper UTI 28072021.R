#######################################################################################
# This is the list of scenarios investigate din the inference paper                   #
# Here we list a set of possible data generating models (dgm) where the true model is #
# indeed one of the family of fractional polynomials:                                 #
#######################################################################################

#For optimise UTI

library(splines)
library(boot)
library(flux)
library(mfp)

#For UUTI we have six durations (4/6/8/10/12/14)

# Scenario 1: Linear relation between logodds and duration:
#formula1.4<-function(x) {
#  alpha<-logit(0.6)
#  beta<-(logit(0.85)-logit(0.6))/4
#  logodds<-alpha+beta*(x-1)
#  return(exp(logodds)/(1+exp(logodds)))
#}
#plot(formula1.4(1:6), type = "l", xlim = c(1,6), ylim=c(0,1))

# Scenario 2: Quadratic +linear relation between logodds and duration:
#formula2.4<-function(x) {
#  alpha<-logit(0.55)
#  beta<-(logit(0.80)-logit(0.55))/4
#  gamma<-(logit(0.85)-logit(0.80))/1
#  logodds<-alpha+beta*(x-1)+gamma*(x-1)^2
#  return(exp(logodds)/(1+exp(logodds)))
#}
#plot(formula2.4(1:6), type = "l", xlim = c(1,6), ylim=c(0,1))

# Scenario 3: Quadratic relation between logodds and duration:
#formula3.4<-function(x) {
#  alpha<-logit(0.6)
#  beta<-(logit(0.85)-logit(0.6))/4^2
#  logodds<-alpha+beta*(x-1)^2
#  return(exp(logodds)/(1+exp(logodds)))
#}
#plot(formula3.4(1:6), type = "l", xlim = c(1,6), ylim=c(0,1))

# Scenario 4: constant 
#formula4.4<-function(x) {
#  alpha<-0.85+0*x
#  return(alpha)
#}
#plot(formula4.4(1:6), type = "l", xlim = c(1,6), ylim=c(0,1))

# Scenario 5: Logarithmic relation between logodds and duration:
#formula5.4<-function(x) {
#  alpha<-logit(0.6)
#  beta<-exp((logit(0.85)-logit(0.6))/4)
#  logodds<-alpha+beta*(log(x))
#  return(exp(logodds)/(1+exp(logodds)))
#}
#plot(formula5.4(1:6), type = "l", xlim = c(1,6), ylim=c(0,1))

# Scenario 6: Square rooted relation between logodds and duration:
#formula6.4<-function(x) {
#  alpha<-logit(0.55)
#  beta<-(logit(0.85)-logit(0.55))/sqrt(4)
#  logodds<-alpha+beta*sqrt(x-1)
#  return(exp(logodds)/(1+exp(logodds)))
#}
#plot(formula6.4(1:6), type = "l", xlim = c(1,6), ylim=c(0,1))

# Scenario 7: Cubic relation between logodds and duration:
#formula7.4<-function(x) {
#  alpha<-logit(0.65)
#  beta<-(logit(0.85)-logit(0.45))/4^3
#  logodds<-alpha+beta*(x-1)^3
#  return(exp(logodds)/(1+exp(logodds)))
#}
#plot(formula7.4(1:6), type = "l", xlim = c(1,6), ylim=c(0,1))

# Scenario 8: Cubic +quadratic relation between logodds and duration:
#formula8.4<-function(x) {
#  alpha<-logit(0.7)
#  beta<-(logit(0.75)-logit(0.7))/4^2
#  gamma<-(logit(0.75)-logit(0.75))/4^3
#  logodds<-alpha+beta*(x-1)^2+gamma*(x-1)^3
#  return(exp(logodds)/(1+exp(logodds)))
#}
#plot(formula8.4(1:6), type = "l", xlim = c(1,6), ylim=c(0,1))

####10/04/2021 Just use these scenarios for now and re-label them 1 to 8

# Scenario 9: Logistic growth model, early growth
formula1.4<-function(x) {
  return(0.05+0.80/(1+exp(-x*2+12)))
}

# Scenario 10: Logistic growth model, growth in the middle
formula2.4<-function(x) {
  return(0.05+0.80/(1+exp(-x*2+18)))
}

# Scenario 11: Gompertz model, b<1:
formula3.4<-function(x) {
  return(0.8*exp(-1*exp(-0.5*(x-6))))
}

# Scenario 12: Gompertz model, b=1:
formula4.4<-function(x) {
  return(0.8*exp(-1*exp(-1*(x-4))))
}

# Scenario 13: Gompertz model, b>1: 
formula5.4<-function(x) {
  return(0.8*exp(-1*exp(-2*(x-3))))
}

# Scenario 14: Quadratic Relation between cure proportion and duration; curvature>0 
formula6.4<-function(x) {
  alpha<-0.6
  beta<-((0.85)-(0.6))/(11)^2
  return(alpha+beta*(x-4)^2)
}

# Scenario 15: Quadratic Relation between cure proportion and duration; curvature<0
formula7.4<-function(x) {
  alpha<-0.6
  gamma<-(0.85-0.6)/(-1*(13)^2)
  beta<-(-2*gamma*(13))
  return(alpha+gamma*(x-4)^2+beta*(x-4))
}

# Scenario 16: linear spline
formula8.4<-function(x) {
  alpha1<-0.65
  alpha2<-0.75
  alpha3<-0.8
  beta1<-(alpha2-alpha1)/3
  beta2<-(alpha3-alpha2)/3
  beta3<-0.001
  p<-rep(0,length(x))
  for (i in 1: length(x)) {
    if (x[i]<4) {
      p[i]<-alpha1+beta1*(x[i]-4)
    } else if (x[i]<6) {
      p[i]<-alpha2+beta2*(x[i]-6)
    } else {
      p[i]<-alpha3+beta3*(x[i]-12)
    }
  }
  return(p)
}

#################################################################################
# This is set of possible data generating models (dgm) used in the first paper: #
#################################################################################

# Scenario 1: Logistic growth model, early growth
#formula1.1<-function(x) {
#  return(0.05+0.80/(1+exp(-x*2+5)))
#}

# Scenario 2: Logistic growth model, growth in the middle
#formula2.1<-function(x) {
#  return(0.05+0.80/(1+exp(-x*2+6)))
#}

# Scenario 3: Logistic growth model, late growth 
#formula3.1<-function(x) {
#  return(0.05+0.80/(1+exp(-x*2+8)))
#}

# Scenario 4: Gompertz model, b<1:
#formula4.1<-function(x) {
#  return(0.8*exp(-1*exp(-0.5*(x-3.5))))
#}

# Scenario 5: Gompertz model, b=1:
#formula5.1<-function(x) {
#  return(0.8*exp(-1*exp(-2*(x))))
#}

# Scenario 6: Gompertz model, b>1: 
#formula6.1<-function(x) {
#  return(0.8*exp(-1*exp(-2*(x))))
#}

# Scenario 7: Linear relation between logodds and duration:
#formula7.1<-function(x) {
#  alpha<-logit(0.7)
#  beta<-(logit(0.95)-logit(0.7))/10
#  logodds<-alpha+beta*(x-10)
#  return(exp(logodds)/(1+exp(logodds)))
#}

# Scenario 8: Quadratic Relation between cure proportion and duration; curvature>0 
#formula8.1<-function(x) {
#  alpha<-0.7
#  beta<-((0.95)-(0.8))/(10)^2
#  return(alpha+beta*(x-10)^2)
#}

# Scenario 9: Quadratic Relation between cure proportion and duration; curvature<0
#formula9.1<-function(x) {
#  alpha<-0.7
#  gamma<-(0.95-0.7)/(-1*(10)^2)
#  beta<-(-2*gamma*(10))
#  return(alpha+gamma*(x-10)^2+beta*(x-10))
#}

# Scenario 10: linear spline
#formula10.1<-function(x) {
#  alpha1<-0.4
#  alpha2<-0.70
#  alpha3<-0.84
#  beta1<-(alpha2-alpha1)/2
#  beta2<-(alpha3-alpha2)/3
#  beta3<-0.001
#  p<-rep(0,length(x))
#  for (i in 1: length(x)) {
#    if (x[i]<12) {
#      p[i]<-alpha1+beta1*(x[i]-1)
#    } else if (x[i]<15) {
#      p[i]<-alpha2+beta2*(x[i]-2)
#    } else {
#      p[i]<-alpha3+beta3*(x[i]-3)
#    }
#  }
#  return(p)
#}


#######################################################################################
# Here we list a set of possible data generating models (dgm) where the true model is #
# indeed one of the family of fractional polynomials:                                 #
#######################################################################################

# Scenario 1: Linear relation between logodds and duration:
#formula1.2<-formula7.1
#plot(formula1.2(1:21), type = "l", xlim = c(10,21), ylim=c(0,1))

# Scenario 2: Quadratic +linear relation between logodds and duration:
#formula2.2<-function(x) {
#  alpha<-logit(0.55)
#  beta<-(logit(0.80)-logit(0.55))/4
#  gamma<-(logit(0.85)-logit(0.80))/1
#  logodds<-alpha+beta*(x-1)+gamma*(x-1)^2
#  return(exp(logodds)/(1+exp(logodds)))
#}
#plot(formula2.2(1:21), type = "l", xlim = c(10,21), ylim=c(0,1))

# Scenario 3: Quadratic relation between logodds and duration:
#formula3.2<-function(x) {
#  alpha<-logit(0.6)
#  beta<-(logit(0.85)-logit(0.6))/4^2
#  logodds<-alpha+beta*(x-1)^2
#  return(exp(logodds)/(1+exp(logodds)))
#}
#plot(formula3.2(1:21), type = "l", xlim = c(10,21), ylim=c(0,1))

# Scenario 4: constant 
#formula4.2<-function(x) {
#  alpha<-0.85+0*x
#  return(alpha)
#}
#plot(formula4.2(1:21), type = "l", xlim = c(10,21), ylim=c(0,1))

# Scenario 5: Logarithmic relation between logodds and duration:
#formula5.2<-function(x) {
#  alpha<-logit(0.6)
#  beta<-exp((logit(0.85)-logit(0.6))/4)
#  logodds<-alpha+beta*(log(x))
#  return(exp(logodds)/(1+exp(logodds)))
#}
#plot(formula5.2(1:21), type = "l", xlim = c(10,21), ylim=c(0,1))

# Scenario 6: Square rooted relation between logodds and duration:
#formula6.2<-function(x) {
#  alpha<-logit(0.55)
#  beta<-(logit(0.85)-logit(0.55))/sqrt(4)
#  logodds<-alpha+beta*sqrt(x-1)
#  return(exp(logodds)/(1+exp(logodds)))
#}
#plot(formula6.2(1:21), type = "l", xlim = c(10,21), ylim=c(0,1))

# Scenario 7: Cubic relation between logodds and duration:
#formula7.2<-function(x) {
#  alpha<-logit(0.65)
#  beta<-(logit(0.85)-logit(0.45))/4^3
#  logodds<-alpha+beta*(x-1)^3
#  return(exp(logodds)/(1+exp(logodds)))
#}
#plot(formula7.2(1:21), type = "l", xlim = c(10,21), ylim=c(0,1))

# Scenario 8: Cubic +quadratic relation between logodds and duration:
#formula8.2<-function(x) {
#  alpha<-logit(0.7)
#  beta<-(logit(0.75)-logit(0.7))/4^2
#  gamma<-(logit(0.75)-logit(0.75))/4^3
#  logodds<-alpha+beta*(x-1)^2+gamma*(x-1)^3
#  return(exp(logodds)/(1+exp(logodds)))
#}
#plot(formula8.2(1:21), type = "l", xlim = c(10,21), ylim=c(0,1))

#######################################################################################
# Here we list a set of possible data generating models (dgm) where the true model is #
# NOT one of the family of fractional polynomials:                                    #
#######################################################################################


# Scenario 1: Logistic growth model, early growth
#formula1.3<-function(x) {
#  return(0.05+0.80/(1+exp(-x*2+5)))
#}

# Scenario 2: Logistic growth model, growth in the middle
#formula2.3<-function(x) {
#  return(0.05+0.80/(1+exp(-x*2+6)))
#}

# Scenario 3: Gompertz model, b<1:
#formula3.3<-function(x) {
#  return(0.8*exp(-1*exp(-0.5*(x-3.5))))
#}

# Scenario 4: Gompertz model, b=1:
#formula4.3<-function(x) {
#  return(0.8*exp(-1*exp(-1*(x-2))))
#}

# Scenario 5: Gompertz model, b>1: 
#formula5.3<-function(x) {
#  return(0.8*exp(-1*exp(-2*(x))))
#}

# Scenario 6: Quadratic Relation between cure proportion and duration; curvature>0 
#formula6.3<-function(x) {
#  alpha<-0.6
#  beta<-((0.85)-(0.7))/(4)^2
#  return(alpha+beta*(x-1)^2)
#}

# Scenario 7: Quadratic Relation between cure proportion and duration; curvature<0
#formula7.3<-function(x) {
#  alpha<-0.6
#  gamma<-(0.85-0.6)/(-1*(4)^2)
#  beta<-(-2*gamma*(4))
#  return(alpha+gamma*(x-1)^2+beta*(x-1))
#}

# Scenario 8: linear spline
#formula8.3<-function(x) {
#  alpha1<-0.4
#  alpha2<-0.70
#  alpha3<-0.84
#  beta1<-(alpha2-alpha1)/3
#  beta2<-(alpha3-alpha2)/3
#  beta3<-0.001
#  p<-rep(0,length(x))
#  for (i in 1: length(x)) {
#    if (x[i]<11) {
#      p[i]<-alpha1+beta1*(x[i]-1)
#    } else if (x[i]<14) {
#      p[i]<-alpha2+beta2*(x[i]-2)
#    } else {
#      p[i]<-alpha3+beta3*(x[i]-3)
#    }
#  }
#  return(p)
#}

##################################################################
# Plot scenarios and check that they are not inferior at 17 days #
##################################################################

# library(boot)
# 
# is.noninf<-matrix(NA,8,2)
# 
# for (i in 2:3) {
#   for (j in 1:8) {
#     plot(get(paste("formula",j,".", i, sep=""))(1:21), type = "l", xlim = c(10,21), ylim=c(0,1), ylab = paste("formula", j,".",i,sep=""))
#     is.noninf[j,i-1]<-(get(paste("formula",j,".", i, sep=""))(20)-get(paste("formula",j,".", i, sep=""))(17)<0.095)
#   }
# }
# apply(is.noninf,2,mean)


