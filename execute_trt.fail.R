x <- rep(0, 10)

p <- 0 # proportion of herd receiving prophylaxis
herd.size <- 1
initial.infection <- 0
vector.pop <- 5000

x[1] <- herd.size * p 
x[2] <- herd.size * (1 -p) - initial.infection
x[4] <- initial.infection 
x[8] <- vector.pop5/12

a <- 0.15/4
b <- 0.48
beta <- a*b

parms <- c(mu.host =      0,   #birth
           mu.vec  =      0.1, #birth
           lambda.host =  0, #death
           lambda.vec =   0.1, #dearh
           
           beta.r = 0,         # infection rate
           beta.s = 0.1,         # infection rate
           epsilon = 0, # treatment failure rate
           
           gamma = 0,     #recovery
           
           delta.r = 0,  #  treatment failure leads to breakthrough
           delta.s = 0,  #  treatment failure leads to breakthrough
           tau = 0, # treatment administered rate
           
           phi = 0, # emergence of resistance
           rho = 0, #resusceptibility
           
           a1 = a,
           c1 = 0.025,
           T1 = 1/4
                      
)

times <- seq(0, 10000, 1)


out <- as.data.frame(ode(x, times, trt.fail.amr, parms))

tail(out)

par(mfrow=c(1,2))

plot(out[,3] ~ out[,1],
     type = "l", col = "blue",
     ylim = c(0,55))
lines(out[,5] ~ out[,1], col = "red")

plot(log(out[,9]) ~ out[,1],
     type = "l", col = "blue",
     ylim = c(0,20))
lines(log(out[,10]) ~ out[,1], col = "red")

