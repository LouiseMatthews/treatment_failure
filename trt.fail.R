
# Load packages ----------------------------------------------------------------

library(deSolve)


# Model function ---------------------------------------------------------------

trt.fail.amr <- function(times, x, parms){
    
  Sp  = x[1]  # Susceptible host with prophylaxis
  Sn  = x[2]  # Susceptible host with no prophylaxis
  Ir  = x[3]  # Infected host with resistant strain
  Is  = x[4]  # Infected host with sensitive strain
  Irt = x[5]  # Infected host with resistant strain w/ treatment
  Ist = x[6]  # Infected host with sensitive strain w/ treatment
  R   = x[7]  # Recovered host
  Sv  = x[8]  # Susceptible vector
  Ivs = x[9]  # Infected vector with sensitive strain
  Ivr = x[10] # Infected vector with resistant strain
  
  
  with(as.list(parms),
       {
         # Calculate population totals
         Nh = Sp + Sn + Ir + Is + Irt + Ist + R   # Number of hosts
         Nv = Sv + Ivs + Ivr                      # Number of vectors 
         
         # Differential equations
          
         dSp.dt  <- (mu.host * Sp - lambda.host * Sp) -
                    (beta.r * (Ivr/Nh) * Sp - (1 - epsilon) * beta.s * (Ivs/Nh) * Sp) +
                    (rho * p * R)
         
         dSn.dt  <- (mu.host * Sn - lambda.host * Sn) -
                    (beta.r * (Ivr/Nh) * Sn - beta.s * (Ivs/Nh) * Sn) +
                    (rho * (1 - p) * R)
         
         dIr.dt  <- ((beta.r * (Ivr/Nh) * Sp) + (beta.r * (Ivr/Nh) * Sn)) -
                    (tau * Ir)
         
         dIs.dt  <- (((1 - epsilon) * beta.s * (Ivs/Nh) * Sp) + (beta.s * (Ivs/Nh) * Sn)) -
                    (tau * Is)
         
         dIrt.dt <- (tau * Ir) - (delta.r * Irt) - (gamma * Irt)
         
         dIst.dt <- (tau * Is) - (delta.s * Ist) - (gamma * Ist) - (phi * Ist)
         
         dR.dt   <- (delta.r * Irt) + (delta.s * Ist) - (rho*p * R) - (rho*(1 - p) * R)
         
         
         dSv.dt  <- (mu.vec * Nv - lambda.vec * Sv)  - 
                    exp(-lambda.vec*T1) * (c1 * a1 * (((Ir + Irt)/Nh) *Sv)) -
                    exp(-lambda.vec*T1) * (c1 * a1 * (((Is + Ist)/Nh) * Sv))
           
         
         dIvs.dt <- exp(-lambda.vec*T1) * (c1 * a1 * (((Is + Ist)/Nh) * Sv)) - (lambda.vec * Ivs)
         
         dIvr.dt <- exp(-lambda.vec*T1) * (c1 * a1 * (((Ir + Irt)/Nh) * Sv)) - (lambda.vec * Ivr)
         
         
         # Define outputs
         dx <- c(dSp.dt, dSn.dt, dIr.dt, dIs.dt, dIrt.dt, dIst.dt, dR.dt, 
                 dSv.dt, dIvs.dt, dIvr.dt)
         list(dx)
       })
  
}

