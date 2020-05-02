
# Load packages ----------------------------------------------------------------

library(deSolve)

# Load function ----------------------------------------------------------------

source("trt.fail.R")



# Run function -----------------------------------------------------------------


x <- rep(0, 10) # set up vector of initial conditions

p <- 0 # proportion of herd receiving prophylaxis

herd.size <- 50           # initial number of cattle
initial.infection <- 0    # initial number of drug sensitive infections
vector.pop <- 5000        # initial population of vectors

# fill initial vector with relevant data
Sp_init  = herd.size * p  # Susceptible host with prophylaxis
Sn_init  = herd.size * (1 -p) - initial.infection  # Susceptible host with no prophylaxis
Ir_init  = 0  # Infected host with resistant strain
Is_init  = initial.infection   # Infected host with sensitive strain
Irt_init = 0  # Infected host with resistant strain w/ treatment
Ist_init = 0  # Infected host with sensitive strain w/ treatment
R_init  = 0  # Recovered host
Sv_init  = vector.pop  # Susceptible vector
Ivs_init = 0  # Infected vector with sensitive strain
Ivr_init = 0
x <- c(Sp = Sp_init, Sn = Sn_init, Ir = Ir_init, Is = Is_init, 
                 Irt = Irt_init, Ist = Ist_init, R = R_init, Sv = Sv_init,
                 Ivs = Ivs_init, Ivr = Ivr_init)


# calculate betas
a <- 0.15/4
b <- 0.46
beta <- a*b

parms <- c(mu.host =      0,   #birth
           mu.vec  =      0, #birth
           lambda.host =  0, #death
           lambda.vec =   0, #dearh
           
           beta.r = 0,         # infection rate
           beta.s = beta,         # infection rate
           epsilon = 0, # treatment failure rate
           
           gamma = 0,     #recovery
           
           delta.r = 0,  #  treatment failure leads to breakthrough
           delta.s = 0,  #  treatment failure leads to breakthrough
           tau = 0, # treatment administered rate
           
           phi = 0, # emergence of resistance
           rho = 0, #resusceptibility
           
           a1 = a,
           c1 = 0.025,
           T1 = 20 # incubation period (days)
                      
)

times <- seq(0, 1000, 1)


out <- as.data.frame(ode(x, times, trt.fail.amr, parms))

tail(out)

par(mfrow=c(1,2))

plot(out$Sn ~ out$time,
     type = "l", col = "blue",
     ylim = c(0,55))
lines(out$Is ~ out$time, col = "red")

plot(out$Sv ~ out$time,
     type = "l", col = "blue",
     ylim = c(0,10000))
lines(log(out$Ivs) ~ out$time, col = "red")

