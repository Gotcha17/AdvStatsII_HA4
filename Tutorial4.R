
# Problem 1 ---------------------------------------------------------------

# Part a)
# define a function for the CRLB of the quantities
CRLB = function(theta, n, c = 0){
  t1 = theta^2/n
  t2 = 1/(n*theta^2)
  t3 = c^2/(n*theta^2)*exp(-2*c/theta)
  return(list(t1, t2, t3))
  # this will return the CRLB's as a list, which is more handy in the plot loop
}

# Part b)
# define a function for the three estimators
estim = function(x, c=0){
  t1 = mean(x)
  t2 = (length(x)-1)/sum(x)
  t3 = mean(as.numeric(x>c))
  return(c(t1, t2, t3))
}

# Part c)
# define global variables
set.seed(666)            # sets the seed for drawing random numbers
MC    = 1000             # number of replications
n     = 1:100            # number of draws/observations in a sample
theta = c(0.1, 1, 2, 5)  # true values of theta
c     = c(1, 10)         # values for the constant c in the third estimator

for (i in 1:3){
# generate one pdf for each estimator (each gets a different Name)
  path=file.path(getwd(),paste("Problem_1_",i,".pdf",sep=""))
  pdf(file = path)
  if (i<=2){
    # for the first two estimators set the figure grid on 2x2
    par(mfrow = c(2,2))
    for (j in theta){
      VAR = c()       # empty vector
      for (k in n){
        # compute the variances of the estimators and store them
        VAR[k]=var(replicate(MC,estim(x=rexp(k,rate=1/j))[i]))
      }
      # plot the variance on the sample size
      # also add appropiate labels and meaningful title
      plot(n,VAR,type="l",col="blue",xlab="n",ylab="variance", 
           ylim=c(0,max(VAR)+0.01),main=bquote(theta==.(j)),lwd=2)
      # add the CRLB to the plot
      lines(n,CRLB(theta=j,n=n)[[i]],col="green",lwd=2)
      # add a meaningfull legend
      legend("topright",legend=c(paste("estimator_",i),"CRLB"), fill=c("blue","green"))
      box()
    }
  }else{
    # for the third estimator set the figure grid on 4x2
    # the rest is the same as before (only one more loop for c)
    par(mfrow = c(4,2))
    for (j in theta){
      for (l in c){
        VAR = c()
        for (k in n){
          VAR[k]=var(replicate(MC,estim(x=rexp(k,rate=1/j),c=l)[i]))
        }
        plot(n,VAR,type="l",col="blue",xlab="n",ylab="variance", 
             ylim=c(0,max(VAR)+0.01),lwd=2, 
             main=substitute(paste(theta," = ",u,", ",c," = ",v), 
                             list(u=j,v=l)))
        lines(n,CRLB(theta=j,n=n,c=l)[[i]],col="green",lwd=2)
        legend("topright",legend=c(paste("estimator_",i),"CRLB"), 
               fill=c("blue","green"))
        box()
      }
    }
  }
  dev.off()
}

# d) Discussion of the results (see gures on the page 4 to 6):

# The CRLB is valid for all estimators.  Only for the second estimator, given a sample of only
# one observation, the variance is always zero and therefore below the CRLB, nevertheless the
# CRLB is only valid for unbiased estimators which is not the case here (why?).  The variance
# of the third estimator in case of small theta
# and large c is also zero which is smaller than the
# CRLB (which is virtually zero itself).
# The  rst estimator reaches the CRLB for every theta
# and n.  The second one for larger n and all theta.
# The last estimator reaches the CRLB for larger n and parameter combinations that
# guarantee that an observation full lls the condition x elenbt fi (c,infinity).


# Problem 2 ---------------------------------------------------------------

# Part a)

# define the log-likelihood
logilog = function(theta, x){
  ll = -length(x)*theta+log(theta)*sum(x)-sum(log(factorial(x)))
  return(ll)
}

# Part b)

# generate a poisson distributed sample of your choice
x = rpois(n = 10000, lambda = 5)

# Part c)

# calculate numerically the MLE (by opimizing the log-likelihood function)
theta_ml_emp = optimize(f=logilog,lower=0.001,upper=10000,x=x,maximum=T)
# and by the analytically derived form (sample mean)
theta_ml_theo = mean(x)

# Part d)

# state the computation time for the numerical method
system.time(optimize(f=logilog,lower=0.001,upper=10000,x=x,maximum=T))

# and for the simple calculation of the mean
system.time(mean(x))


# Problem 3 ---------------------------------------------------------------

# Part a)

# again define the log-likelihood
logilog2 = function(mu, x){
  ll = length(x)*log(1/sqrt(2*pi))-1/2*sum(x^2)+mu*sum(x)-length(x)/2*mu^2
  return(ll)
}

# Part b)

# define global variables
n = c(10, 100, 1000, 10000)
mu_theo = seq(from = -5, to = 5, length.out = 1000)
# simulation and plot
path = file.path(getwd(), paste("Problem_3.pdf"))
pdf(file = path)
par(mfrow = c(2,2))
start = proc.time() # save computation time till the loop in variable "start"
for (i in n){
  mu_emp = c()
  k      = 1
  for (j in mu_theo){
    # generate sample
    x = rnorm(i, mean = j, sd = 1)
    # calculate numerically the MLE and store it
    mu_emp[k] = optimize(f=logilog2,lower=-100,upper=100,x=x,maximum=TRUE)[[1]]
    k = k+1
  }
  plot(mu_theo, mu_emp, pch = '.', col = "red", main = bquote(n == .(i)))
  lines(mu_theo, mu_theo, col = "green")
}
dev.off()

proc.time()-start # computation time of the loop (c)

# Discussion of the results (see figures on the page 9):
# With increasing n the theta_hat_ML approach the true value of theta.
