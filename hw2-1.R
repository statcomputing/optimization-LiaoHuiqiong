#1(b)
#graph the log-likelihood function
x <- c(1.77, -0.23, 2.76, 3.80, 3.47, 56.75, -1.34, 4.24, -2.44,
       3.29, 3.71, -2.40, 4.53, -0.07, -1.05, -13.87, -2.53, -1.75)
th <- seq(-50,50,0.5)
n <- length(x)
ll <- c()
for (i in 1:length(th)) {
  ll[i] <- -n*log(pi)-sum(log(1+(th[i]-x)^2))
}
plot(th,ll,type = "l", ylab = "ll", xlab = "theta", main = "Cauchy Log Likelihood ")

#find MLE for theta using the Newton-Raphson method
cauchy.ll.1st.deriv <- function(theta){
  z <- -2*(sum((theta-x)/(1+(theta-x)^2)))
  return(z)
}
cauchy.ll.2nd.deriv <- function(theta){
  y <- -2*(sum((1-(theta-x)^2)/((1+(theta-x)^2)^2)))
  return(y)
} 
staring.p <- c(-11,-1, 0, 1.5, 4, 4.7, 7, 8, 38)
sample_mean <- mean(staring.p)
print(sample_mean)

NewtonR <- function(theta0){
  theta <- array()
  theta[1] <- theta0
  i <- 1
  difference <- 10
  while (abs(difference)>10^(-10) & i<200){
    theta[i+1] <- theta[i] - cauchy.ll.1st.deriv(theta[i]) / cauchy.ll.2nd.deriv(theta[i])
    difference <- abs(theta[i+1] - theta[i])
    i <- i + 1 
  }
  return(theta[i])
}
NewtonR(-11)
NewtonR(-1)
NewtonR(0)
NewtonR(1.5)
NewtonR(4)
NewtonR(7)
NewtonR(8)
NewtonR(38)
NewtonR(sample_mean)

# When starting at -11,7,8,and 38, the function cannot find the optimial point
# Sample mean(5.68889) also is not a good starting point

#1(c)
#alpha = 1, 0.64 or 0.25
#alpha <- -1 / cauchy.ll.2nd.deriv(theta)
alpha1 <- 1 
fixP1 <- function(theta0){
  theta <- array()
  theta[1] <- theta0
  i <- 1
  difference <- 10
  while (abs(difference)>10^(-10) & i<200) {
    theta[i+1] <- theta[i] + alpha1 * cauchy.ll.1st.deriv(theta[i])
    difference <- abs(theta[i+1] - theta[i])
    i <- i + 1
  }
  return(theta[i])
}

fixP1(-11)
fixP1(-1)
fixP1(0)
fixP1(1.5)
fixP1(4)
fixP1(4.7)
fixP1(7)
fixP1(8)
fixP1(38)

alpha2 <- 0.64
fixP2 <- function(theta0){
  theta <- array()
  theta[1] <- theta0
  i <- 1
  difference <- 10
  while (abs(difference)>10^(-10) & i<200) {
    theta[i+1] <- theta[i] + alpha2 * cauchy.ll.1st.deriv(theta[i])
    difference <- abs(theta[i+1] - theta[i])
    i <- i + 1
  }
  return(theta[i])
}

fixP2(-11)
fixP2(-1)
fixP2(0)
fixP2(1.5)
fixP2(4)
fixP2(4.7)
fixP2(7)
fixP2(8)
fixP2(38)

alpha3 <- 0.25
fixP3 <- function(theta0){
  theta <- array()
  theta[1] <- theta0
  i <- 1
  difference <- 10
  while (abs(difference)>10^(-10) & i<200) {
    theta[i+1] <- theta[i] + alpha3 * cauchy.ll.1st.deriv(theta[i])
    difference <- abs(theta[i+1] - theta[i])
    i <- i + 1
  }
  return(theta[i])
}

fixP3(-11)
fixP3(-1)
fixP3(0)
fixP3(1.5)
fixP3(4)
fixP3(4.7)
fixP3(7)
fixP3(8)
fixP3(38)

#1(d)use fisher score to find MLE for theta

fisher.score <- n / 2
FisherFun <- function(theta0){
  theta <- array()
  theta[1] <- theta0
  i <- 1
  difference <- 10
  while (abs(difference)>10^(-10) & i<200){
    theta[i+1] <- theta[i] + cauchy.ll.1st.deriv(theta[i]) / (fisher.score)
    difference <- abs(theta[i+1] - theta[i])
    i <- i + 1
  }
  return(theta[i])
}
FisherFun(-11) 
FisherFun(-1)
FisherFun(0)
FisherFun(1.5)
FisherFun(4)
FisherFun(4.7)
FisherFun(7)
FisherFun(8)
FisherFun(38)
#when FisherFun starts at -11, -1, and 0, it will return an optimal value at -0.5914735
#when FisherFun starts at 1.5, 4, 4.7, 7, and 8, it will return an optimal value at 3.021345
#when starts at 38, it will return at 3.037648

#Refine: plug -0.5914735, 3.021345, and 3.037648 into NewtonR()
NewtonR(-0.5914735)
NewtonR(3.021345)
NewtonR(3.037648)

#1(e)