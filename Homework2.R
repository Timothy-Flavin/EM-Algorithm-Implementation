#Stochastic Modeling Homework 2. Normal pdf


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
xpts = seq(0,6,0.1)
plot(type = 'l', x=xpts, y=GaussMixtureSeq(xpts, c(0.2,0.2,0.2,0.4), c(1,2,3,4), c(.1,.1,.1,.1) ) )

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
    # draw from dstribution i
    sample[point] <- rnorm(n=1,mean=us[i],sd=sqrt(vars[i]))
  }
  return(sample)
}

EM <- function(pts, pis, us, vars, epsilon){
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
        w[d] <- dnorm(pts[i],mean=temp.us[d],sd=sqrt(temp.vars[d])) * temp.pis[d]
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
    
    print("--------------------------------------------------")
    print("pis: ")
    print(old.pis)
    print(temp.pis)
    print("uss: ")
    print(old.us)
    print(temp.us)
    print("vars: ")
    print(old.vars)
    print(temp.vars)
    
    if(plotProgress){
      plot(type = 'l', x=xpts, y=GaussMixtureSeq(xpts, temp.pis, temp.us, temp.vars ) )
    }
    
    max.dif = max(  c( max(abs(temp.pis-old.pis), max(abs(temp.us-old.us)), max(abs(temp.vars-old.vars))) )  )
    print(max.dif)
  }
}

EM(sample, c(0.2,0.2,0.2,0.4), c(0,2,2.5,3), c(.15,.15,.15,.15), epsilon=0.001 )


c()
c(1)
append(c(), 2)
append(c(1), 2)

sample <- GaussMixtureSample(10000, c(0.2,0.2,0.2,0.4), c(1,2,3,4), c(.1,.1,.1,.1))

hist(sample, breaks=100)

plot(type = 'l', x=xpts, y=GaussMixtureSeq(xpts, c(0.20207,0.20449,0.22747,0.36597), c(0.9945625,2.0318276,3.0861390,4.0596749), c(0.08859049,0.07284075,0.08331663,0.07280922) ) )
sample <- GaussMixtureSample(100000, c(0.20207,0.20449,0.22747,0.36597), c(0.9945625,2.0318276,3.0861390,4.0596749), c(0.08859049,0.07284075,0.08331663,0.07280922))
hist(sample, breaks=50)
