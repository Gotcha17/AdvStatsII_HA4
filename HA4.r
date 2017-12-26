# Define variables
mu_theo = seq(from = -5, to = 5, length.out = 100)
sigma_theo = seq(from = 0.1, to = 25, length.out = 100)
n = seq(from = 100, to = 1000, length.out = 10)

##### a) ##################################################################################################

# Log-likelihood function for the normal dirstribution
# mu and sigma are inputed through the "par" vector, this is needed due to the use of the optim() fct.
loglike = function(par, x){
  mu = par[1]
  sigma = par[2]
  ll = -length(x)/2*log(2*pi*sigma^2)-1/(2*sigma^2)*sum(x^2)+mu/(sigma^2)*sum(x)-length(x)/(2*sigma^2)*mu^2
  return(ll)
}

##### b) ##################################################################################################

# Euclidean distance function for two points/vectors
euclidean_dist = function(p, q){
  d = sqrt(sum((p-q)^2))
  return(d)
}

##### c) ##################################################################################################

# Simulation
mu_emp = array(dim = c(100,100,10))                           # 3-dimenional R-object is initialized
sigma_emp = array(dim = c(100,100,10))                        # 3-dimenional R-object is initialized
d = 1                                                         # initializing index for level dimension
for (i in n){                                                 # for-loop over all n
  r = 1                                                       # initializing index for row dimension
  for (j in mu_theo){                                         # for-loop over all theoretical mus
    c = 1                                                     # initializing index for column dimension
    for (l in sigma_theo){                                    # for-loop over all theoretical sigmas
      # generate sample
      x = rnorm(i, mean = j, sd = l)                          # generating n nor. dist. r.v.
      # calculate numerically the MLE and store it
      # starting params: 0 for mu and 1 for sigma, maximize(fnscale=-1) fct. loglike, with given data x 
      optim = optim(par=c(0,1), fn=loglike, x=x, control=list(fnscale=-1)) 
      mu_emp[r,c,d] = optim$par[1]                            # stores empirical value for mu in 3-D obj.
      sigma_emp[r,c,d] = optim$par[2]                         # stores empirical value for mu in 3-D obj.
      c = c+1                                                 # increasing index for column dimension
    }
    r = r+1                                                   # increasing index for row dimension
  }
  d = d+1                                                     # increasing index for level dimension
}

##### d) ##################################################################################################

# Calculating the euclidean distance for every possible combination of mu and sigma for empirical and theo.
eu_dist_array = array(dim = c(100,100,10))                    # 3-dimenional R-object is initialized
for (d in seq(from = 1, to = 10)){                            # for-loop over all n
  for (r in seq(from = 1, to = 100)){                         # for-loop over all mus
    for (c in seq(from = 1, to = 100)){                       # for-loop over all sigmas
      # Eucl. dist. is calculated with formula from part b) and results are stored in a 3-D R-object
      eu_dist_array[r,c,d] = euclidean_dist(c(mu_emp[r,c,d],sigma_emp[r,c,d]),c(mu_theo[r],sigma_theo[c]))
    }
  }
}

##### e) ##################################################################################################

# In order to see whether there is a dependence on the starting value, we examine the graph of the
# log-likelihood function as a function of mu and sigma for a given random sample
# as shown in the plots, the log-likelihood function of the normal distribution is (globally) concave 
# and has a unique maximum. Therefore the used optim() function (using the Nelder and Mead method) 
# achieves this maximum independent of given starting values. Therefore there is no dependence of 
# euclidean distance and the starting values.
# However this is not the case for any maximum likelihood (c.f. normal mixture) where starting values 
# may have a great impact on results obtained.
# Also one might state, that there could be an influence of the maximum iteration steps and the
# optimization method used.

#Plot
set.seed(123)
x = rnorm(100000, mean=0, sd=1)

# compute log-likelihood as function of mu
Mu = seq(-1,1, length.out=100)
LL1 = rep(0,length(Mu))
for (i in 1:length(LL1)) {
  LL1[i] = loglike(c(Mu[i],1), x=x)
}
optim = optim(par=c(0,1), fn=loglike, x=x, control=list(fnscale=-1)) 
plot(Mu, LL1, type="l", col="blue", lwd=2)
abline(v=optim$par[1], col="red", lwd=2)

# compute log-likelihood as function of sigma
Sigma = seq(0.5,2, length.out=100)
LL2 = rep(0,length(Sigma))
for (i in 1:length(LL2)) {
  LL2[i] = loglike(c(0,Sigma[i]), x=x)
}

plot(Sigma, LL2, type="l", col="blue", lwd=2)
abline(v=optim$par[2], col="red", lwd=2)

###########################################################################################################
