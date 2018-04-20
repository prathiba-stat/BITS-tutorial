library(runjags) 
#setwd('C:/enter working directory here and remove the # sign')
source("plot_SSD.R")
y <- c(78, 83, 86, 118, 72, 116, 129, 177, 152, 121,
       93, 93, 100, 112, 122, 176, 156, 159, 147, 135) 
yre <- matrix(y, 4, 5, byrow = TRUE) #reformat the data for plotting
#reformat into 4 rows and 5 columns meaning 4 phases and 5 observations/phase
plot_SSD(yre) #plot the data 
T <- length(y) #total number of observations
P <- 4  #no. of phases
Tp <- c(5, 5, 5, 5)     #length of each phase
Tt <- c(0, cumsum(Tp))        #cumulative sum tells us the time-point
beta <- c() #empty vector to store phase means as starting values
for (l in 1:P){
  beta <- c(beta, mean(y[(Tt[l] + 1):Tt[l + 1]]))
}

#Model is defined:
BITS.model1.4 <- "model {
    for (l in 1:P){
      yhat[(Tt[l] + 1)] <- beta[l]
      y[(Tt[l] + 1)] ~ dnorm(yhat[(Tt[l] + 1)], tau) 
        for (i in (Tt[l] + 2):Tt[l + 1]){
          yhat[i] <- beta[l]
          y[i] ~ dnorm(yhat[i] + rho * (y[i - 1] - yhat[i - 1]), tau)
        }
        beta[l] ~ dnorm(mu[l], 0.01)
        mu[l] ~ dnorm(40, .05)
    }
   sigma ~ dunif(0.1, 5)
   tau <- pow(sigma, -2)
   rho ~ dunif(-1, 1)
}"
#end of model definition

#Begin running the model with the data
results <- autorun.jags(
  model = BITS.model1.4, 
  data = list(y = y, Tt = Tt, P = P), 
  monitor = c("beta", "sigma", "rho", "es"), 
  n.chains = 4,
  startsample = 30000,
  inits = function() { 
    list(
      beta = apply(as.matrix(beta), 1, function(x)rnorm(1, x, 1)),
      sigma = runif(1, 0.1, 5),
      rho = runif(1, -1, 1)
    )
  }, 
  method = "rjparallel"
)


# combine all chains into a single chain for convenience
results$draws <- combine.mcmc(results$mcmc)
results.agg <- data.frame(cbind(results$HPD, results$summary$statistics))
write.csv(results.agg, "ABAB-results.csv")


