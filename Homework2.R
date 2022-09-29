# Stochastic Modeling Homework 2. Normal pdf

# ################## FILE SUMMARY ############################################
#                                                                            #
# This file contains function definitions first and then code that executes  #
# them near the bottom. The functions are GaussMixtureSeq which returns the  #
# pdf values for a sequence of x values given some mixture parameters.       #
# GaussMixtureSample returns a sample from a mixture with 'n' points in it.  #
# OldEM is a neutered version of the EM algorith that for some reason is way #
# more stable. EM() is the em algorithm using the update rules from this     #
# homework. The code that executes parts b,c and d is at the bottom.         #
#                                                                            #
##############################################################################



# returns pdf value at each point in pts using a mixture
# model defined by sum(pis[i] * N(pts[i],us[i],vars[i])
# Where N(pt, u, var) returns the value of the pdf defined
# by a normal distribution with mean u and variance var 
GaussMixtureSeq <- function (pts, pis, us, vars){
  result <- numeric(length(pts))
  for(i in 1:length(pts)){
    for(p in 1:length(pis)){
      result[i] <- result[i] + pis[p] * dnorm(pts[i],mean=us[p],sd=sqrt(vars[p]))
    }
  }
  return(result)
}

GaussMixtureSample <- function(n, pis, us, vars){
  # pick a random distribution, i, from the pi's
  sample = numeric(length=n)
  for(point in 1:n){
    z = runif(n=1,min=0,max=1)
    i = 1
    while(z>0){
      z<-z-pis[i]
      i<-i+1
    }
    i=i-1
    # draw from dstribution[i]
    sample[point] <- rnorm(n=1,mean=us[i],sd=sqrt(vars[i]))
  }
  return(sample)
}

# This is an old EM algorithm that I used a long time ago that does a few things
# sort of wrong, but for some reason it is MUCH more stable than the "correct" 
# algorithm so I am going to leave it here to ask you about it before class
OldEM <- function(pts, pis, us, vars, epsilon){
  groups <- list(length=length(pis))
  
  temp.pis <- pis
  temp.us <- us
  temp.vars <- vars
  
  max.dif = 1.0
  while(max.dif > epsilon && min(temp.pis) > 0.01){
    plotProgress = TRUE
    for(i in 1:length(pis)){
      groups[[i]] <- c(1)
    }
    for(i in 1:length(pts)){
      w <- numeric(length(pis))
      for(d in 1:length(pis)){
        w[d] <- dnorm(pts[i],mean=temp.us[d],sd=sqrt(temp.vars[d])) #* temp.pis[d]
      }
      groups[[which.max(w)]] <- append(groups[[which.max(w)]], pts[i])
    }
    
    old.pis <- temp.pis
    old.us <- temp.us
    old.vars <- temp.vars
    
    for(i in 1:length(pis)){
      groups[[i]] <- groups[[i]][-1]
      temp.pis[i] <- 1.0*length(groups[[i]])/length(pts)
      temp.us[i] <- mean(groups[[i]])
      temp.vars[i] <- var(groups[[i]])
    }
    
    if(plotProgress){
      graph.xs<-seq(floor(min(pts)),floor(max(pts))+1,by=0.01)
      plot(type = 'l', x=graph.xs, y=GaussMixtureSeq(graph.xs, temp.pis, temp.us, temp.vars ) )
    }
    
    max.dif = max(  c( max(abs(temp.pis-old.pis), max(abs(temp.us-old.us)), max(abs(temp.vars-old.vars))) )  )
  }
  params = list("pis"=temp.pis, "us"=temp.us,"vars"=temp.vars)
}


# Problem II part c.
# pts: vector of points sampled from the mixture density
# pis: also called alphas, the vector of proportions
# us: the vector of means for the distributions. Should be same length as pis
# vars: vector of variances. should be same length as pis
EM <- function(pts, pis, us, vars, epsilon){
  responsibilities = matrix(0, nrow=length(pis), ncol=length(pts), byrow=TRUE)
  num.iter <- 0
  temp.pis <- pis
  temp.us <- us
  temp.vars <- vars
  temp.probs <- c(1:length(pis)) # Just setting aside some memory to hold p_i's
  max.dif = 1
  # When none of the parameters change much, or one of the proportions or 
  # variances is nearly zero, the algorithm exits. The variance or proportions
  # go to zero when the EM algorithm decides that there are less densities than
  # were originally estimated
  while(max.dif > epsilon && min(temp.pis) > epsilon && min(temp.vars)>epsilon/100 && num.iter<1000){
    num.iter<- num.iter+1
    old.pis <- temp.pis
    old.us <- temp.us
    old.vars <- temp.vars
    
    # Graphs the current estimate at each iteration to play a nice little animation
    graph.xs<-seq(floor(min(pts)),floor(max(pts))+1,by=0.01)
    plot(type = 'l', x=graph.xs, y=GaussMixtureSeq(graph.xs, temp.pis, temp.us, temp.vars ) )
    
    # rix.sums will hold the vector of sums of r_i(k) * x_k to keep from 
    # backtracking  
    rix.sums <- numeric(length(pis))
    
    # ri.sums will hold the vector of sums r_i(k) for each i to save for later
    ri.sums <- numeric(length(pis))
    for(k in 1:length(pts)){
      for(i in 1:length(pis)){
        temp.probs[i] <- dnorm(pts[k],old.us[i],old.vars[i]) * old.pis[i]
      }
      rk.denom <- sum(temp.probs)
      for(i in 1:length(pis)){
        # Calculating the responsibilities matrix for look-up later because 
        # memory is cheap these days
        responsibilities[i,k] <- temp.probs[i] * 1.0 / rk.denom
        
        # might as well record these sums ill re-use while im at it
        ri.sums[i] <- ri.sums[i] + responsibilities[i,k]
        rix.sums[i] <- rix.sums[i] + responsibilities[i,k]*pts[k]
      }
    }
    
    # Each of these lines is taken from the update rules defined in Module 8b
    # and also found in part (II a.) of this homework
    # Time to update the parameters:
    for(i in 1:length(pis)){
      temp.pis[i] <- 1.0/length(pts) * ri.sums[i]
      temp.us[i] <- 1.0/ri.sums[i] * rix.sums[i]
      temp.vars[i] <- 1.0/ri.sums[i] * (responsibilities[i,] %*% (pts-temp.us[i])^2)
    }
    
    # This finds the parameter with the biggest change since the last update
    max.dif = max(  c( max(abs(temp.pis-old.pis), max(abs(temp.us-old.us)), max(abs(temp.vars-old.vars))) )  )
  }
  params = list("pis"=temp.pis, "us"=temp.us,"vars"=temp.vars)
  return(params)
}

# Problem II part b.
set.seed(42)
alpha <- c(0.1,0.2,0.3,0.4)
mu <- c(2,3,5,7)
vars <- c(0.08,0.1,0.2,0.1)
xpts = seq(0,10,0.01)

plot(type = 'l', x=xpts, y=GaussMixtureSeq(xpts, alpha, mu, vars ) )
sample <- GaussMixtureSample(10000, alpha, mu, vars)
hist(sample, breaks=50)

vars<-vars*2
plot(type = 'l', x=xpts, y=GaussMixtureSeq(xpts, alpha, mu, vars ) )


# Problem II part c/d: Show the EM Algorithm working
# Example 1
set.seed(42)
alpha <- c(0.1,0.2,0.3,0.4)
mu <- c(1,3,5,8)
vars <- c(0.08,0.1,0.1,0.1)
xpts = seq(0,10,0.01)

plot(type = 'l', x=xpts, y=GaussMixtureSeq(xpts, alpha, mu, vars ) )
sample <- GaussMixtureSample(10000, alpha, mu, vars)
hist(sample, breaks=50)
params <- EM(sample, c(0.11,0.3,0.2,0.39), c(1.1,2.8,5.1,8), c(.4,.4,.4,0.4), epsilon=0.0001 )
params

# Example 2
set.seed(42)
alpha <- c(0.7,0.1,0.05,0.15)
mu <- c(1,3,5,8)
vars <- c(0.08,0.08,0.1,0.1)
xpts = seq(0,10,0.01)

plot(type = 'l', x=xpts, y=GaussMixtureSeq(xpts, alpha, mu, vars ) )
sample <- GaussMixtureSample(10000, alpha, mu, vars)
hist(sample, breaks=50)
params <- EM(sample, c(0.11,0.3,0.2,0.39), c(1.1,2.8,5.1,8), c(.4,.4,.4,0.4), epsilon=0.0001 )
params

# This is an example where the em algroithm does some weird crap where it 
# decides that some of the clusters actually don't exist and instead there are
# only 2. The OldEM that is more like kmeans works still
set.seed(42)
sample <- GaussMixtureSample(10000, c(0.2,0.2,0.2,0.4), c(1,3,4,5), c(.1,.05,.05,.1))
hist(sample, breaks=100)
params <- EM(sample, c(0.21,0.21,0.21,0.37), c(0.9,3.1,4.1,6), c(.8,.5,.5,1.0), epsilon=0.001 )
print(params)
params <- OldEM(sample, c(0.21,0.21,0.21,0.37), c(0.9,3.1,4.1,6), c(.8,.5,.5,1.0), epsilon=0.001 )
print(params)
params <- EM(sample, params[["pis"]], params[["us"]], params[["vars"]], epsilon=0.001 )
print(params)

# Moving the mean of the big guy over by only 1 makes it work again with this seed
# sometimes it still doesn't work on other seeds so this case must be close to bad
set.seed(42)
sample <- GaussMixtureSample(10000, c(0.2,0.2,0.2,0.4), c(1,3,4,6), c(.1,.05,.05,.1))
hist(sample, breaks=100)
params <- EM(sample, c(0.21,0.21,0.21,0.37), c(0.9,3.1,4.1,6), c(.8,.5,.6,1.0), epsilon=0.0001 )
print(params)
params <- OldEM(sample, c(0.21,0.21,0.21,0.37), c(0.9,3.1,4.1,6), c(.8,.5,.6,1.0), epsilon=0.0001 )
print(params)
params <- EM(sample, params[["pis"]], params[["us"]], params[["vars"]], epsilon=0.0001 )
print(params)
