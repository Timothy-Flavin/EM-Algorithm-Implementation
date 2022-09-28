#Stochastic Modeling Homework 2. Normal pdf


# returns pdf value at each point in pts using a mixture
# model defined by sum(pis[i] * N(pts[i],us[i],vars[i])
# Where N(pt, u, var) returns the value of the pdf defined
# by a normal distribution with mean u and variance var 
GaussMixture <- function (pts, pis, us, vars){
  result <- numeric(length(pts))
  for(i in 1:length(pts)){
    for(p in 1:length(pis)){
      result[i] <- result[i] + pis[p] * dnorm(pts[i],mean=us[p],sd=sqrt(vars[p]))
    }
  }
  return(result)
}
xpts = seq(0,10,0.1)
plot(type = 'l', x=xpts, y=GaussMixture(xpts, c(0.1,0.3,0.2,0.4), c(1,2,3,5), c(.1,.1,.1,.3) ) )
