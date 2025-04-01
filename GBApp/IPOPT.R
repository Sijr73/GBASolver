# GBA optimization for general N (formulation on scaled fluxes w)

library(ipoptr)


# Functions ################################################################################

# mu ############## (growth rate)
mu <- function(f) as.numeric(M[p,r]*f[r]/(tau(ci(f))%*%f)) 

# negative mu (for minimization)
neg_mu <- function(f) -as.numeric(M[p,r]*f[r]/(tau(ci(f))%*%f)) 
# gt ############## (growth timje = 1/mu)
#gt <- function(w) as.numeric(tau(ci(w))%*%w)/(N[a,r]*w[r])

# fluxes
v <- function(f) as.vector(mu(f))*rho*f

# protein concentrations
prot <- function(f) tau(ci(f))*v(f)

# internal concentrations "i" of metabolites "m" and total protein "a"
ci <- function(f) rho*M%*%f
E <- function(f) rho*dtau(ci(f))%*%M
# scaled concentrations
b <- function(f) M%*%f
# Define the objective function and its derivative
neg_mu <- function(f) -as.numeric(M[p, r] * f[r] / (tau(ci(f)) %*% f))
neg_dmu <- function(f) -((mu(f)^2) / b(f)[p]) * (M[p, ] / mu(f) - t(f) %*% (rho * dtau(ci(f)) %*% M) - tau(ci(f)))
Gj <- function(f,j) M[p,j] + mu(f)*(- tau(ci(f))[j] - t(f)%*%(E(f)[,j])
                                    
                                    + c(t(f)%*%E(f)%*%f)*sM[j])  
f <- function(fy) c(1/sM[1] - sM[-1]%*%fy,fy)

# f1 as a function of fy #######################################################
f1 <- function(fy) (1 - sM[-1]%*%fy)/sM[1]


# Define the equality constraint function and its derivative

g <- function(f) sM %*% f - 1
dg <- function(w) sM

# Define the inequality constraint function and its derivative
h <- function(f) c(ci(f) , rho * tau(ci(f)) * f)
dh	<- function(f) rho*rbind(M,rho*(dtau(ci(f))%*%M*matrix(rep(f,r),nrow=r)) + 
                              diag(r)*tau(ci(f)) )
# Define the IPOPT callback function for the objective, constraints, and derivatives
eval_f <- function(f) neg_mu(f)
eval_grad_f <- function(f) neg_dmu(f)
eval_g <- function(f) as.numeric(c(g(f), h(f)))
eval_jac_g <- function(f) t(rbind(dg(f), dh(f)))

#initial points
source('/app/f0_shiny.R',local = TRUE)

# Define IPOPT options

# Solve the optimization problem using IPOPT
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
  #source('w0.R')
  
  rho <- rho_cond[cond]
  
  x  <- x_cond[,cond]
  
  # Upper bounds 
  upper_f <- rep(10,r)
  
  # lower bounds 
  lower_f <- rep(-10,r)
  
  
  constraint_lb <- c(0, rep(0, length(eval_g(f0)) - 1))
  constraint_ub <- c(0, rep(Inf, length(eval_g(f0)) - 1))
  eval_jac_g_structure <- rep(list(c(1:r)),(length(eval_g(f0))))
  
  opts <- list(
    "tol" = 1.0e-5,
    "max_iter" = 1000,
    "linear_solver" = "ma57",
    #"file_print_level"=12,
    #"output_file"="banana.out",
    "derivative_test"="first-order",
    #"jacobian_approximation" = "finite-difference-values",
    "gradient_approximation"="exact"
     #"warm_start_init_point" = "yes"
  )
  
  
  # Optimization
  
  #   measuring the total optimization time
  st <- system.time({
    
    
    res <- ipoptr(x0 = f0,
                  eval_f = eval_f,
                  eval_grad_f = eval_grad_f,
                  eval_g = eval_g,
                  eval_jac_g = eval_jac_g,
                  eval_jac_g_structure = eval_jac_g_structure,
                  constraint_lb = constraint_lb,
                  constraint_ub = constraint_ub,
                  opts = opts)
  })
  
  
  # solution
  # w_opt[cond,] <- res$par
  f0 <- res$solution
  #f_opt[cond,] <- res$solution
  f_opt[cond,] <- f0
  
  # optimal mu
  mu_opt[cond] <- mu(f0)
  

  # convergence (codes: "-1" means optimization problem, "5" means optimization stop because maxeval, "4" means
  # optimization stopped because xtol_rel or xtol_abs (above) was reached)
  
  #conv[cond] <- res$convergence
  conv[cond] <- res$status
  
  # optimization time
  otime[cond] <- signif(st[[1]], digits= 4)
  
  # number of iterations
  #iter[cond] <- res$iter
  iter[cond] <- res$iterations
  
  
  print(paste("optimization: ",cond, "/", n_conditions,", optimization time: ",otime[cond]," s, convergence: ",
              conv[cond],", growth rate: ",signif(mu_opt[cond], digits= 3), sep=""))
  
  
  # lambda
  lambda[cond] <- c(t(f0)%*%E(f0)%*%f0)*(mu(f0)^2)/b(f0)[p]
  
  # Growth adaptation coefficient with respect to the density rho
  A_rho[cond]  <- -lambda[cond]/mu(f0)
  
  # if(res$status==0){
  #    w0 <- res$solution
  #  }
  
  for (j in 1:r) dmu_opt[cond,j] <- Gj(f0,j)
  incProgress(1/n_conditions, detail = paste("Condition", cond))
  
  # Pause for 0.1 seconds to simulate a long computation.
  Sys.sleep(0.001)
  
}
  # next initial value
  
})

# produces c_opt
c_opt <- matrix(rep(0,p*n_conditions),ncol=p)
for (cond in 1:n_conditions) c_opt[cond,] <- ci(f_opt[cond,])

# mean optimization time
mean_time <- signif(mean(otime), digits= 3)

dmu_opt
