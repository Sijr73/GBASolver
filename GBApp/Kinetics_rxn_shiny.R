# Turnover time function for different kinetic rate laws, and their partial
# derivatives with respect to c
#
# NOTE: external concentrations x are treated as global variables, so they are 
# not inputs in the functions


# irreversible MM ##############################################################
iMM <- function(j,c,xc) as.numeric(prod(1 + KS[,j]/xc)/kcatf[j]) 
    

# derivative with respect to c
diMM <- function(j,c,xc){
  
  ditauj <- rep(0,p)
  
  for (i2 in 1:p) {
      
      y <- i2 + nx # position of i in terms of all reactants xc
      
      ditauj[i2] <- -as.numeric((KS[y,j]/(c[i2]^2))*prod(1 + 
                                                             KS[-y,j]/xc[-y])/kcatf[j]) 
      
    }
    
  return(ditauj)
  
}

# derivative with respect to x

xdiMM <- function(j,c,xc){
  
  dntauj <- rep(0,nx)
  
  for (n in 1:nx) {
      
      dntauj[n] <- -as.numeric((KS[n,j]/(x[n]^2))*prod(1 + 
                                                           KS[-n,j]/xc[-n])/kcatf[j]) 
      
    }
  
  return(dntauj)
  
}

# irreversible MM + inhibition (here only one inhibitor per reaction) ##########

# (it's useful to first define the matrix of reciprocal KI, with zero for KI=0)

rKI <- 1/KI
rKI[rKI == Inf] <- 0

iMMi <- function(j,c,xc) as.numeric(prod(1 + xc*rKI[,j])*prod(1 + KS[,j]/xc)/kcatf[j]) 
    
# derivative with respect to ci

diMMi <- function(j,c,xc){
  
  ditauj <- rep(0,p)
  
  for (i2 in 1:p) {
      
      y <- i2 + nx # position of i in terms of all reactants
      
      ditauj[i2] <- as.numeric(rKI[y,j]*prod(1 + KS[,j]/xc) - 
                                   prod(1 + xc*rKI[,j])*(KS[y,j]/(c[i2]^2))*prod(1 + KS[-y,j]/xc[-y]))/kcatf[j] 
      
  }
  
  return(ditauj)
  
}

# irreversible MM + activation (only one activator  per reaction) #############

iMMa <- function(j,c,xc) as.numeric(prod(1 + KA[,j]/xc)*prod(1 + KS[,j]/xc)/kcatf[j]) 
    
# derivative with respect to c

diMMa <- function(j,c,xc){
  
  ditauj <- rep(0,p)
  
    for (i2 in 1:p) {
      
      y <- i2 + nx # position of i in terms of all reactants
      
      ditauj[i2] <- -as.numeric( (KA[y,j]/(c[i2]^2))*prod(1 + KS[,j]/xc) + 
                                     (KS[y,j]/(c[i2]^2))*prod(1 + KA[,j]/xc)*prod(1 + KS[-y,j]/xc[-y]) )/kcatf[j] 
      
    }
    
  return(ditauj)
  
}


# derivative with respect to x #################################################

xdiMMa <- function(j,c,xc){
  
  dntauj <- rep(0,nx)
  
  for (n in 1:nx) {
    
    dntauj[n] <- -as.numeric( (KA[n,j]/(x[n]^2))*prod(1 + KS[,j]/xc) + 
                                (KS[n,j]/(x[n]^2))*prod(1 + KA[,j]/xc)*prod(1 + KS[-n,j]/xc[-n]) )/kcatf[j] 
    
  }
  
  return(dntauj)
  
}

# irreversible MM + inhibition + activation ####################################################

iMMia <- function(j,c,xc)  as.numeric(prod(1 + xc*rKI[,j])*prod(1 + KA[,j]/xc)*prod(1 + KS[,j]/xc)/kcat[j]) 
    
# derivative with respect to c (assuming only one inhibitor and activator per reaction) 

diMMia <- function(j,c,xc){
  
  ditauj <- rep(0,p)
  
    for (i2 in 1:p) {
      
      y <- i2 + nx # position of i in terms of all reactants
      
      ditauj[i2] <- as.numeric( rKI[y,j]*prod(1 + KA[,j]/xc)*prod(1 + KS[,j]/xc) + 
                                    
                                    prod(1 + xc*rKI[,j])*(-KA[y,j]/(c[i2]^2))*prod(1 + KS[,j]/xc) + 
                                    
                                    prod(1 + xc*rKI[,j])*prod(1 + KA[,j]/xc)*(-KS[y,j]/(c[i2]^2))*prod(1 + KS[-y,j]/xc[-y]) )/kcat[j] 
      
    }
  
  return(ditauj)
  
}



# Reversible MM  ######################################################################################

# first separate Km matrices for substrates and for products
KS <- K*(Mtotal<0)
KP <- K*(Mtotal>0)

rMM <- function(j,c,xc) as.numeric( 1 / ( kcatf[j]/prod(1 + KS[,j]/xc) - kcatb[j]/prod(1 + KP[,j]/xc)  ) )
    

# derivative with respect to c #################################################

drMM <- function(j,c,xc){
  
  ditauj <- rep(0,p)
  
      for (i in 1:p) {
      
      y <- i + nx # position of i in terms of all reactants xc
      
      # this will be multiplied my -(tau(c)^2) after loops
      ditauj[i] <- as.numeric( (kcatf[j]/prod(1 + KS[-y,j]/xc[-y]))*KS[y,j]/((c[i] + KS[y,j])^2) 
                                 
                                 
                                 -(kcatb[j]/prod(1 + KP[-y,j]/xc[-y]))*KP[y,j]/((c[i] + KP[y,j])^2) ) 
      
    }
    
    ditauj <-  -ditauj*as.numeric( 1 / ( kcatf[j]/prod(1 + KS[,j]/xc) - kcatb[j]/prod(1 + KP[,j]/xc)  )^2 )
  
  return(ditauj)
  
}

# derivative with respect to x #################################################

xdrMM <- function(j,c,xc){
  
  ditauj <- rep(0,nx)
  
      for (n in 1:nx) {
      
      ditauj[n] <- as.numeric( (kcatf[j]/prod(1 + KS[-n,j]/xc[-n]))*KS[n,j]/((x[n] + KS[n,j])^2) 
                                 
                                 
                                 -(kcatb[j]/prod(1 + KP[-n,j]/xc[-n]))*KP[n,j]/((x[n] + KP[n,j])^2) )
      
    }
    
    ditauj<-  -ditauj*as.numeric( 1 / ( kcatf[j]/prod(1 + KS[,j]/xc) - kcatb[j]/prod(1 + KP[,j]/xc)  )^2 )
    
  
  return(ditauj)
  
}

# Finds simplest kinetic functions #############################################

kinetics <- rep(0,r)
for (j in 1:r) {
  
  
  if (kcatb[j] == 0) {
    # irreversible kinetics
    
    if (sum(KI[,j]) == 0 & sum(KA[,j]) == 0 ) kinetics[j] <- paste("iMM(",j,",c,xc)",sep="")
    
    if (sum(KI[,j]) == 0 & sum(KA[,j]) > 0 ) kinetics[j] <- paste("iMMa(",j,",c,xc)",sep="")
    
    if (sum(KI[,j]) > 0 & sum(KA[,j]) == 0 ) kinetics[j] <- paste("iMMi(",j,",c,xc)",sep="")
    
    if (sum(KI[,j]) > 0 & sum(KA[,j]) > 0 ) kinetics[j] <- paste("iMMia(",j,",c,xc)",sep="")
    
  } else {
    
    kinetics[j] <- paste("rMM(",j,",c,xc)",sep="")
    
  }
  
}

# derivative functions
dkinetics <- paste("d",kinetics,sep="")

# derivatives with respect to external concentrations (epsilon)
xdkinetics <- paste("x",dkinetics,sep="")

# Global kinetic function ######################################################

tau <- function(c) {
  
  xc <- c(x,c)
  
  tauj <- rep(0,r)
  
  for (j in 1:r) tauj[j] <- eval(parse(text=noquote(kinetics[j])))
  
  tauj
  
}

# Global kinetic derivatives function ##########################################

dtau <- function(c) {
  
  xc <- c(x,c)
  
  dtauj <- matrix(rep(0,r*p),ncol=p)
  
  for (j in 1:r) dtauj[j,] <- eval(parse(text=noquote(dkinetics[j])))
  
  dtauj
  
}

# Global epsilon function ######################################################

epsilon  <- function(c) {
  
  xc <- c(x,c)
  
  dtauj <- matrix(rep(0,r*nx),ncol=nx)
  
  for (j in 1:r) dtauj[j,] <- eval(parse(text=noquote(xdkinetics[j])))
  
  dtauj
  
}



