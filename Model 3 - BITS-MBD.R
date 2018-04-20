library(runjags)
#setwd('C:/enter working directory here and remove the # sign')
source('plot_SSD_mbd.R')
#Model is defined:
BITS.model1 <- "model {
  yhat[1] <- beta[1, 1]
  yhat[(Tb + 1)] <- beta[2, 1]
  for (i in 2:Tb) {
    yhat[i] <- beta[1, 1]
    y[i] ~ dnorm(yhat[i] + rho * (y[i - 1] - yhat[i - 1]), tau)
  }
  for (i in (Tb + 2):T) {
    yhat[i] <- beta[2, 1]
    y[i] ~ dnorm(yhat[i] + rho * (y[i - 1] - yhat[i - 1]), tau)
  }
  y[1] ~ dnorm(yhat[1], tau) 
  y[(Tb + 1)] ~ dnorm(yhat[(Tb + 1)], tau) 
  es <- (beta[2, 1] - beta[1, 1])/sigma
  for (i in 1:P){
    beta[i, 1] ~ dnorm(mu[i], 0.01)
    mu[i] ~ dnorm(40, .05)
  }
  sigma ~ dunif(0.1, 5)
  tau <- pow(sigma, -2)
  rho ~ dunif(-1, 1)
}"
############# end of model definition##########################
#data is a dataframe with the first, second, and third columns 
#indicating the subject, the values, and the phase
data <- data.frame(cbind(c(rep(1, 10), rep(2, 16), rep(3, 20)),
             c(78, 83, 86, 118, 72, 
              116, 129, 177, 152, 121,
              70, 91, 110, 101, 96, 95, 134, 94, 
              106, 138, 123, 155, 135, 134, 164, 115,
              66, 113, 121, 106, 135, 72, 70, 105, 131, 96, 
              122, 96, 136, 100, 137, 114, 100, 124, 109, 111),
              c(rep(c(1, 2), each = 5),
              rep(c(1, 2), each = 8),
              rep(c(1, 2), each = 10)))) 
colnames(data) <- c("subject", "DV", "phase")
data <- cbind(data, c(seq(1:10), seq(1:16), seq(1:20))) #reformat data
colnames(data) <- c("case", "outcome", "treatment", "time") #column names
data$treatment <- data$treatment - 1 
plot_SSD_mbd(data) #plot data
P <- 2
#change nSubject to reflect the number of subjects in your data
#change nbase to reflect the baselengths in your data
#change nTime to reflect the lengths in your data
nSubject <- 3         #There are 3 subjects in this dataset
nbase <- c(5, 8, 10)  #Baseline length for each subject
nTime <- c(10, 16, 20) #total length for each subject
pcalc.agg <- matrix(c(1:nSubject, rep(NA, nSubject)), 
                    nSubject, 2) #empty matrix to store pcalcs
colnames(pcalc.agg) <- c("subject", "p.calc")
results.agg <- data.frame()
#begin for loop with one iteration per subject
for (i in 1:nSubject){
  y <- data[data$subject==i,2]
  T <- nTime[i]
  Tb <- nbase[i]
  beta1 <- mean(y[1:Tb])
  beta2 <- mean(y[(Tb + 1):T])
  results <- autorun.jags(
    model = BITS.model1, 
    data = list(y = y, T = T, P = P, Tb = Tb), 
    monitor = c("beta", "sigma", "rho", "es"), 
    n.chains = 3,
    startsample = 30000,
    inits = function() { 
      list(
        beta = rbind(rnorm(1, beta1, 1), rnorm(1, beta2, 1)),
        sigma = runif(1, 0.1, 5),
        rho = runif(1, -1, 1)
      )
    }, 
    method = "rjparallel"
  )
  results$draws <- combine.mcmc(results$mcmc)
  results.tab <- summary(results)
  pcalc.agg[i,2] <- length(intersect(results$draws[,"beta[1,1]"], 
                             results$draws[,"beta[2,1]"]))/length(results$draws[,"beta[1,1]"])
  results.agg <- rbind.data.frame(results.agg, data.frame(cbind("subject" = i, results$HPD, 
                                  results$summary$statistics)))
}
write.csv(results.agg, 'MBD-results.csv')
write.csv(pcalc.agg, 'p-calc-MBD.csv', row.names = FALSE)

