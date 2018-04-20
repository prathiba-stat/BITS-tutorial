# plot a single-subject design

plot_SSD <- function(y) {
  require(ggplot2)
  
  P <- nrow(y)  
  T <- apply(y, 1, function(x) max(which(!is.na(x))))
  
  observations <- rep(NA, sum(T))
  k = 1
  for (i in 1:P) {
    for (j in 1:T[i]) {
      observations[k] = y[i, j]
      k = k + 1
    }
  }
  changepoints <- cumsum(T)
  
  DF <- data.frame(
    phase = rep(seq(P), times = T),
    time = seq(sum(T)),
    value = observations
  )
  
  ggplot(DF, aes(x = time, y = value)) +
    geom_point() +
    geom_line(aes(group = phase)) +
    geom_vline(xintercept = head(changepoints, -1) + 0.5, linetype = "dotted") +
    theme_bw()
}