# Growth optimization

# Functions ####################################################################

# mu ############## (growth rate)
mu <- function(f) as.numeric(M[p,r]*f[r]/(tau(ci(f))%*%f)) 

# negative mu (for minimization)
negative_mu <- function(f) -as.numeric(M[p,r]*f[r]/(tau(ci(f))%*%f)) 

# fluxes
v <- function(f) as.vector(mu(f))*rho*f

# protein concentrations
prot <- function(f) tau(ci(f))*v(f)

# internal concentrations "i" of metabolites "m" and proteome "p"
ci <- function(f) rho*M%*%f

# biomass fractions
b <- function(f) M%*%f

# Optimization #################################################################

# equality constraints (density constraint)

g <- function(f) sM%*%f - 1

# inequality constraints (non negativity of concentrations). Here we assume ####
# mu ~ 1 for a faster calculation of p, instead of using the function p(f) #####

h <- function(f) c(ci(f), rho*tau(ci(f))*f)

# Indirect elasticities ########################################################

E <- function(f) rho*dtau(ci(f))%*%M

# Gradient of negative_mu ######################################################

negative_dmu <- function(f) -((mu(f)^2)/b(f)[p])*(M[p,]/mu(f) - 
                                      t(f)%*%(rho*dtau(ci(f))%*%M) - tau(ci(f)))

# Jacobian of equality constraint ##############################################

dg	<- function(f) sM

# Jacobian of inequality constraints ###########################################

# rows = c1, c2,...,p1,p2... , columns = d/df1, d/df2, ...
dh	<- function(f) rho*rbind(M,rho*(dtau(ci(f))%*%M*matrix(rep(f,r),nrow=r)) + 
                              diag(r)*tau(ci(f)) )

# Equation of motion for reaction j (Gj = 0 at optimality) #####################
Gj <- function(f,j) M[p,j] + mu(f)*(- tau(ci(f))[j] - t(f)%*%(E(f)[,j])
      
      + c(t(f)%*%E(f)%*%f)*sM[j])  

# f as a function of the independent fy ########################################
f <- function(fy) c(1/sM[1] - sM[-1]%*%fy,fy)

# f1 as a function of fy #######################################################
f1 <- function(fy) (1 - sM[-1]%*%fy)/sM[1]


# Finds initial f0 #############################################################

source('/app/f0_shiny.R',local = TRUE)

# Gamma y ######################################################################

#Gamma_y <- function(fy) (1 - sM[y]/sM[1])*dmudfy(fy)


# dmu/dfy ######################################################################

#dmudfy <- function(fy) 




# starts loop for the optimization at each environmental condition #############
f_opt  <- matrix(rep(0,r*n_conditions),ncol=r)
mu_opt <- rep(0,n_conditions)
otime  <- rep(0,n_conditions)
conv   <- rep(0,n_conditions)
iter   <- rep(0,n_conditions)
lambda <- rep(0,n_conditions)
A_rho  <- rep(0,n_conditions)
dmu_opt  <- matrix(rep(0,r*n_conditions),ncol=r)
withProgress(message = 'Running simulation', value = 0, {
for (cond in 1:n_conditions) {
  
  rho <- rho_cond[cond]
  
  x  <- x_cond[,cond]
  
  # Upper bounds 
  upper_f <- rep(1000,r)
  
  # lower bounds 
  lower_f <- rep(-1000,r)
  
 
  # Optimization
  #solver      <- "SLSQP" # Local solvers: "COBYLA", "LBFGS", "MMA", or "SLSQP"
  
  # measuring the total optimization time
  st <- system.time({
    
    # optimization (by default the function auglag minimizes, so we use 
    # negative_mu instead of mu as the objective function)
    
    res <- auglag(f0, fn = negative_mu, gr = negative_dmu, lower = lower_f, 
                  upper = upper_f, heq = g, heqjac = dg, hin = h,hinjac = dh, 
                  localsolver = c(solver), localtol = 1e-8,
                  control = list(maxeval = 10000))
    
  })
  
  # solution
  f0 <- res$par
  
  # solution
  f_opt[cond,] <- f0
  
  # optimal mu
  mu_opt[cond] <- mu(f0)
  
  # convergence (codes: "-1" means optimization problem, "5" means optimization 
  # stop because maxeval, "4" means optimization stopped because xtol_rel or 
  # xtol_abs (above) was reached)
    conv[cond] <- res$convergence
  
  # optimization time
  otime[cond] <- signif(st[[1]], digits= 4)
  
  # number of iterations
  iter[cond] <- res$iter  
  
  print(paste("optimization: ",cond, "/", n_conditions,", optimization time: ",
              otime[cond]," s, convergence: ", conv[cond],", growth rate: ",
              signif(mu_opt[cond], digits= 3), sep=""))
  
  # lambda
  lambda[cond] <- c(t(f0)%*%E(f0)%*%f0)*(mu(f0)^2)/b(f0)[p]
  
  # Growth adaptation coefficient with respect to the density rho
  A_rho[cond]  <- -lambda[cond]/mu(f0)
  
  # Gj at fy0
  for (j in 1:r) dmu_opt[cond,j] <- Gj(f0,j)
  incProgress(1/n_conditions, detail = paste("Condition", cond))
  
  # Pause for 0.1 seconds to simulate a long computation.
  Sys.sleep(0.001)
  
}

})
# produces c_opt
c_opt <- matrix(rep(0,p*n_conditions),ncol=p)
for (cond in 1:n_conditions) c_opt[cond,] <- ci(f_opt[cond,])

# mean optimization time
mean_time <- signif(mean(otime), digits= 3)

dmu_opt
