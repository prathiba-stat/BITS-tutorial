#text that follows the # symbol are comments and 
#will not be executed
install.packages(runjags) #this command installs the runjags package
library(runjags) #loads runjags into the environment
#setwd("C:/enter working directory here and remove the comment sign")
#y is the dependent variable with 10 observations, 5 per phase
y <- c(78, 83, 86, 118, 72, 116, 129, 177, 152, 121) 
#T is the total length of y, i.e. 10
T <- length(y)
#P is assigned the value 2, i.e. the number of phases
P <- 2
#Baseline has 5 observations. Change 5 to any number accordind to your data
Tb <- 5
#the following two lines compute the means of the two phases
#beta1 is the mean of the first Tb observations and beta2 - the rest
#these betas will be used as starting values in autorunjags

beta1 <- mean(y[1:Tb])
beta2 <- mean(y[(Tb + 1):T])
#Model is defined:
BITS.model1 <- "model {
#yhat for time-point 1 is set to the intercept or the estimated mean of the first phase
#this is because yhat is the expected value of y and the expected value of y is the mean
    yhat[1] <- beta[1, 1]
#similarly the sixth time-point (i.e. the first time-point of phase 2)
#is assigned to the estimated mean of the second phase
    yhat[(Tb + 1)] <- beta[2, 1]
#for loop begins at 2 and runs till Tb (for baseline phase)
     for (i in 2:Tb) {
#the expected value of yhat in the baseline is beta[1, 1]
        yhat[i] <- beta[1, 1]
#y is drawn from a distribution with expected value yhat 
#and autocorrelation rho which is multiplied by the error of the 
#previous time-point (i - 1)
#tau is the precision = 1/(standard deviation) = 1/sigma
        y[i] ~ dnorm(yhat[i] + rho * (y[i - 1] - yhat[i - 1]), tau)
     }
#for loop begins at Tb+1 and runs till T (for treatment phase)
     for (i in (Tb + 2):T) {
        yhat[i] <- beta[2, 1]
        y[i] ~ dnorm(yhat[i] + rho * (y[i - 1] - yhat[i - 1]), tau)
     }
      y[1] ~ dnorm(yhat[1], tau) #equation 1 for baseline
      y[(Tb + 1)] ~ dnorm(yhat[(Tb + 1)], tau) #equation 1 for treatment
      es <- (beta[2, 1] - beta[1, 1])/sigma
#Prior specifications
#For both phases
      for (i in 1:P){
#the intercepts for baseline and treatment phases are drawn 
#from distributions with means mu[1] and mu[2], respectively 
# and precision 0.01
        beta[i, 1] ~ dnorm(mu[i], 0.01)

#because mu is the number of vocal responses in five-minute intervals, 
#the expected value is set at 40 with a standard deviation of 20 (precision = 1/20 = .05)
#change the values inside () to reflect your belief or literature
#You can set the lower limit of mu to 0 by adding I(0, ) after dnorm(40, .05)
        mu[i] ~ dnorm(40, .05)
      }
#standard deviation can uniformly vary between 0.1 to 5. 
#Change values inside the 
#parentheses to reflect your belief or literature
   sigma ~ dunif(0.1, 5)
#tau is sigma^(-2)
   tau <- pow(sigma, -2)
#autocorrelation can vary uniformly between -1 and 1
   rho ~ dunif(-1, 1)
}"
############# end of model definition##########################

#Begin running the model with the data

results <- autorun.jags(
  #autorun.jags automatically runs the model until convergence is indicated
  model = BITS.model1, #the model is BITS.model1 defined above
  data = list(y = y, T = T, P = P, Tb = Tb), #input data are the y vector of observations,
  #the total number of time-points T, the number of baseline observations Tb, and phases P
  monitor = c("beta", "sigma", "rho", "es"), 
  n.chains = 3,#these are the parameters to be monitored i.e., checked for convergence and obtained
  startsample = 30000,
  #burn-in/discard the first 30000 parameter estimates to ensure that 
  #there is no effect of the starting values on the retained estimates
  inits = function() { #initialize/assign starting values
    #change the specs to see how starting values affect the estimates
    #before and after burning-in. Once burned-in they should not.
    list(
      beta = rbind(rnorm(1, beta1, 1), rnorm(1, beta2, 1)),
      #intercept starting values around the phase means
      sigma = runif(1, 0.1, 5),
      #standard deviation can be any value between -1 and 5
      rho = runif(1, -1, 1)
      #autocorrelations between -1 and 1
    )
  }, 
  method = "rjparallel"
  #run the chains in parallel
)


# combine all chains into a single chain for convenience
results$draws <- combine.mcmc(results$mcmc)

results$HPD["beta[1,1]",]
results
head(results$draws)
traceplot(results$mcmc)

sample.results <- run.jags(
  #autorun.jags automatically runs the model until convergence is indicated
  model = BITS.model1, #the model is BITS.model1 defined above
  data = list(y = y, T = T, P = P, Tb = Tb), #input data are the y vector of observations,
  #the total number of time-points T, the number of baseline observations Tb, and phases P
  monitor = c("beta", "sigma", "rho"), 
  n.chains = 3,#these are the parameters to be monitored i.e., checked for convergence and obtained
  burnin = 0,
  sample = 10000,
  adapt = 1,
  #startsample = 30,
  #burn-in/discard the first 30000 parameter estimates to ensure that 
  #there is no effect of the starting values on the retained estimates
  inits = function() { #initialize/assign starting values
    #change the specs to see how starting values affect the estimates
    #before and after burning-in. Once burned-in they should not.
    list(
      beta = rbind(rnorm(1, beta1, 1), rnorm(1, beta2, 1)),
      #intercept starting values around the phase means
      sigma = runif(1, 0.1, 5),
      #standard deviation can be any value between -1 and 5
      rho = runif(1, -1, 1)
      #autocorrelations between -1 and 1
    )
  }, 
  method = "rjparallel"
  #run the chains in parallel
)
library(MCMCpack)
traceplot(results$mcmc)
densityplot.mcmc(results$mcmc)
x11()
par(mfrow=c(2,1))
plot(head(sample.results$mcmc[,2]), trace = FALSE, density = TRUE, smooth = FALSE,
     xlab = " ", main = " ", col = c("black", "black", "black"), 
     lty = c("solid", "dotted", "dashed"))
plot(results$mcmc[,5], trace = TRUE, density = TRUE, smooth = FALSE,
     xlab = " ", main = " ", col = c("black", "black", "black"), 
     lty = c("solid", "dotted", "dashed"))
plot(density(results$mcmc[,5]))

results.tab <- summary(results)
p.calc <- length(intersect(results$draws[,"beta[1,1]"], results$draws[,"beta[2,1]"]))/length(results$draws[,"beta[1,1]"])

##############################
#######EFFECT SIZE#####################
######################################################
lims <- c(min(c(results$draws[,"beta[1,1]"], results$draws[,"beta[2,1]"])),
          max(c(results$draws[,"beta[1,1]"], results$draws[,"beta[2,1]"])))
x11()
par(mar=c(3.5,3.5,2,1),mgp=c(2,0.7,0), mfcol = c(2, 1))
plot(density(results$draws[,"beta[1,1]"]), main = "Intercept Phase 1",
     xlab = expression(beta[11]), ylab =" ", xlim = lims)
abline(v=results$HPD["beta[1,1]", c(1,3)], lty = "dashed")
plot(density(results$draws[,"beta[2,1]"]), main = "Intercept Phase 2",
     xlab = expression(beta[21]), ylab =" ", xlim = lims)
abline(v=results$HPD["beta[2,1]", c(1,3)], lty = "dashed")

plot(density(results$draws[,"es"]), main = "Effect Size",
     xlab = " ", ylab =" ")
abline(v=results.tab["es", c(1, 3)], lty = "dashed")

##############Small ROPE with 1 SD ABOVE 0
#in the plot below the mid-point of comparison for ROPE is 0.5
#the interval ranges from 0.5-0.5 to 0.5 +0.5
#change this value based on the hypothesis for ROPE
plot(density(results$draws[,"rho"]), main = "Autocorrelation",
     xlab = "rho", ylab =" ")
plot(density(results$draws[,"sigma"]), main = "Standard deviation",
     xlab = "SD", ylab =" ")