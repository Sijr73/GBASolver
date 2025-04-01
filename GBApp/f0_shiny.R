# Initial value subproblem: linear optimization to find maximal ribosome flux
# fraction f^r, with a minimal production of each metabolite. The constraints 
# are mass conservation (M*f = b) and surface flux balance (sM*f = 1)

# minimal biomass fractions b

rho <- rho_cond[1]

# average Km per substrate
aKm <- colSums(t(K[-n,]))

aKm[aKm == 0] <- 10

lowerb <- aKm[-p]/1000


# defining parameters
objective.fn <- c(rep(0,r-1),1)
const.mat <- rbind(M[-p,],sM)
const.dir <- c(rep(">=",p-1),"=")
const.rhs <- c(lowerb, 1)

# solving model

lp.solution <- lp("max", objective.fn, const.mat, const.dir, const.rhs)

# solution to the linear optimization

f0 <- lp.solution$solution

# stops if linear optimization didn't find solution

if (sum(f0) == 0)  {

# second method:

# Moore-Penrose inverse of internal M
P <- ginv(M)

# estimated b0 for metabolites
b0 <- 0
b0[1:(p-1)] <- aKm[-p]/sum(aKm[-p])*0.2
b0[p] <- 0.2

# initial v0 estimated assuming no use of null modes
f0 <- P%*%b0

f0
}
# stops if solution is not a valid state for the nonlinear problem due to
# negative growth rate or concentrations

x  <- x_cond[,1]

#if (TRUE %in% c(mu(f0) < 0, prot(f0) < 0, ci(f0) < 0)) stop("initial f0 is not a
#                                   valid state for the nonlinear optimization")
