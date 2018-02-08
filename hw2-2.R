x2 <- c(3.91, 4.85, 2.28, 4.06, 3.70, 4.04, 5.46, 3.53, 2.28, 1.96,
       2.53, 3.88, 2.22, 3.47, 4.82, 2.46, 2.99, 2.54, 0.52)
staring_p <- c(-11,-1, 0, 1.5, 4, 4.7, 7, 8, 38)
#2(a)
theta <- seq(-pi, pi, 0.01)
n2 <- length(x2)
llh <- rep(0,length(theta))
for (i in 1:length(theta)) {
  llh[i] <-sum(log((1 - cos(x2 - theta[i])) / (2 * pi)))
}
plot(theta, llh, type = "l", ylab = "llh", xlab = "theta", main = "Log Likelihood distribution")

#2(b)
mme <- mean(x2)
print(mme)

#2(c)   
llh1stDeriv <- function (theta){
  w <- -(sum(sin(x2 - theta)/(1 - cos(x2 - theta))))
  return(w)
}
llh2ndDeriv <- function(theta){
  k <- -(sum(1 / (1 - cos(x2 - theta))))
  return(k)
}

nr_method <- function(theta0){
  theta <- array()
  theta[1] <- theta0
  i <- 1
  difference <- 10
  while (abs(difference)>10^(-10) & i<200) {
    theta[i+1] <- theta[i] - llh1stDeriv(theta[i]) / llh2ndDeriv(theta[i])
    difference <- theta[i+1] - theta[i]
    i <- i + 1
  }
  return(theta[i])
}

nr_method(mme)
# When the starting point = mean(x2), the function is converged at theta = 3.170715

#2(d)
nr_method(-2.7)
nr_method(2.7)
# When the starting point = -2.7, the function is converged at theta = -2.668857
# When the starting point = 2.7, the function is converged at theta = 2.848415

#2(e)

new.starting.p <- seq(-pi,pi,length.out = 200)
calculation <- sapply(new.starting.p,nr_method)
converge.point <- round(calculation,6)
counting.result <- as.data.frame(table(converge.point))
print(counting.result)
