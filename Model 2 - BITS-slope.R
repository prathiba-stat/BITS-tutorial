#setwd('C:/enter working directory here and remove the # sign')
source('plot_SSD.R')
library(runjags)
y <- c(85, 85, 83, 83, 82, 83, 84, 85, 84, 84, 
       94, 96, 96, 98, 98, 97, 97, 98, 100, 100)
T <- length(y)
P <- 2
Tb <- 10
beta1 <- mean(y[1:Tb])
beta2 <- mean(y[(Tb + 1):T])
#Model is defined:
BITS.model2 <- "model {
  yhat[1] <- beta[1, 1]
  yhat[Tb+1] <- beta[2, 1]
#baseline phase
    for (i in 2:Tb) {
#slope beta[1, 2] is added to the equation
      yhat[i] <- beta[1, 1] + beta[1, 2] * i
      y[i] ~ dnorm(yhat[i] + rho * (y[i - 1] - yhat[i - 1]), tau)
    }
#treatment phase
    for (i in (Tb + 2):T) {
#slope is multiplied by (t-tb) as given in equation 7 such that
#(i-Tb)  = x represents the xth time-point in the intervention phase
      yhat[i] <- beta[2, 1] + beta[2, 2] * (i - Tb)
      y[i] ~ dnorm(yhat[i] + rho * (y[i - 1] - yhat[i - 1]), tau)
    }
    y[1] ~ dnorm(yhat[1], tau)
    y[(Tb + 1)] ~ dnorm(yhat[(Tb + 1)], tau)
#Prior specifications
    for (i in 1:P){
#the slopes for baseline and treatment phases are drawn 
#from distributions with means mu.s[1] and mu.s[2], respectively 
# and precision 0.01
      beta[i, 1] ~ dnorm(mu[i], 0.01)
      beta[i, 2] ~ dnorm(mu.s[i], 0.01)
      mu[i] ~ dnorm(40, 0.5)
      mu.s[i] ~ dnorm(0, 0.1)
    }
    sigma ~ dunif(0.1, 5)
    tau <- pow(sigma, -2)
    rho ~ dunif(-1, 1)
  }"
############# end of model definition##########################

#Begin running the model with the data

results <- autorun.jags(
  model = BITS.model2, 
  data = list(y = y, T = T, P = P, Tb = Tb), 
  monitor = c("beta", "sigma", "rho"), 
  n.chains = 4,
  startsample = 30000,
  inits = function() { 
    list(
      beta = rbind(c(rnorm(1, beta1, 1), 1), 
                   c(rnorm(1, beta2, 1), 1)),
      sigma = runif(1, 0.1, 5),
      rho = runif(1, -1, 1)
    )
  }, 
  method = "rjparallel"
  #run the chains in parallel
)


# combine all chains into a single chain for convenience
results$draws <- combine.mcmc(results$mcmc)

results
results.tab <- summary(results)
p.calc.int <- length(which(results$draws[,"beta[1,1]"]>min(results$draws[,"beta[2,1]"])))/length(results$draws[,"beta[1,1]"])
p.calc.slope <- length(which(results$draws[,"beta[1,2]"]>min(results$draws[,"beta[2,2]"])))/length(results$draws[,"beta[2,2]"])
int.lims <- c(min(c(results$draws[,"beta[1,1]"], results$draws[,"beta[2,1]"])),
          max(c(results$draws[,"beta[1,1]"], results$draws[,"beta[2,1]"])))
slope.lims <- c(min(c(results$draws[,"beta[1,2]"], results$draws[,"beta[2,2]"])),
              max(c(results$draws[,"beta[1,2]"], results$draws[,"beta[2,2]"])))
x11()
par(mar=c(3.5,3.5,2,1),mgp=c(2,0.7,0), mfcol = c(2, 2))
plot(density(results$draws[,"beta[1,1]"]), main = "Intercept Phase 1",
     xlab = expression(beta[11]), ylab =" ", xlim = int.lims)
abline(v=results.tab["beta[1,1]", c(1, 3)], lty = "dashed")
plot(density(results$draws[,"beta[2,1]"]), main = "Intercept Phase 2",
     xlab = expression(beta[21]), ylab =" ", xlim = int.lims)
abline(v=results.tab["beta[2,1]", c(1, 3)], lty = "dashed")

plot(density(results$draws[,"beta[1,2]"]), main = "Slope Phase 1",
     xlab = expression(beta[12]), ylab =" ", xlim = slope.lims)
abline(v=results.tab["beta[1,2]", c(1, 3)], lty = "dashed")
plot(density(results$draws[,"beta[2,2]"]), main = "Slope Phase 2",
     xlab = expression(beta[22]), ylab =" ", xlim = slope.lims)
abline(v=results.tab["beta[2,2]", c(1, 3)], lty = "dashed")

