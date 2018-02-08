#3(a)
beetle <- data.frame(
  days = c(0, 8, 28, 41, 63, 69, 97, 117, 135, 154),
  beetles = c(2, 47, 192, 256, 768, 896, 1120, 896, 1184, 1024))

pop.growth.fun <- function(t, K, r){
  f <- (2 * K) / (2 + (K - 2) * exp(-r * t))
  return(f)
} 

Gauss.N <- nls(beetles~pop.growth.fun(days, K, r), data = beetle, start = list(K = 1300, r = 0.3), trace = TRUE)
summary(Gauss.N)
plot(days, beetles)
lines(days, predict(Gauss.N))

#(b)
days = c(0, 8, 28, 41, 63, 69, 97, 117, 135, 154)
beetles = c(2, 47, 192, 256, 768, 896, 1120, 896, 1184, 1024)
SSE <- function(K, r){
  sse.value <- sum((beetles - (2*K/(2+(K-2)*exp(-r*days))))^2)
  return(sse.value)
}
SSE.matri <- matrix(0,100,100, byrow = TRUE)
for (i in 1:100){
  for (j in 1:100){
    K <- 300+12*i
    r <- 0+0.003*j
    SSE.matri[i, j] <- SSE(K, r)
  }
}
K <- seq(300, 1500, length.out = 100)
r <- seq(0, 0.3, length.out = 100)
contour(K, r, SSE.matri, col = "blue", main = "Contour Plot", xlab = "Population", ylab = 'Growth Rate')




  