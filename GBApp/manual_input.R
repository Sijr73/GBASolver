# Function to handle file input
processmanual <- function(m1,m2,m3,m4,m5,m6) {

  # Position of parameters
  posM          <- m1
  posK          <- m2
  posKI         <- m3
  posKA         <- m4
  poskcat       <- m5
  posconditions <- m6
 
  
  
  
 
  

  # Getting data from sheets
  
  # Mass fraction matrix Mtotal including external reactants
  Mtotal <- as.matrix(posM)
  
  # reaction and reactant names
  reaction <- colnames(Mtotal)
  
  # Michaelis constant matrix K 
 # if (any(posK[which(Mtotal < 0)]  > 0)) {
      
    K <- as.matrix(posK)
    
  #} else K <- 0.1*(Mtotal<0)
  
  # inhibition constant matrix KI
  if (any(posKI > 0)) {
    
    KI <- as.matrix(posKI)
    
  } else KI <- 0*K
  
  # activation constant matrix KA
  if (any(posKA > 0)) {
    
    KA <- as.matrix(posKA)
    
  } else KA <- 0*K
  
  
  # kcat
  kcatf <- as.numeric(poskcat[1,])
  kcatb <- as.numeric(poskcat[2,])
    names(kcatf)=reaction
  names(kcatb)=reaction
  kcat_matrix <- rbind(kcatf, kcatb)
  kcatf=kcat_matrix[1,]
  kcatb=kcat_matrix[2,]
  # Growth condition names
  condition <- colnames(posconditions)
  
  # Mass density rho at each condition
  rho_cond <- as.numeric(posconditions[1,])
  
  # external concentrations at each condition
  x_cond  <-as.matrix(posconditions[-1,])
  if (ncol(x_cond) == 1) {
    x_cond <- t(x_cond)
  }
  conditiontab =rbind(rho_cond,x_cond)
  rownames(conditiontab)=rownames(posconditions)
  rho_cond=conditiontab[1,]
  x_cond=conditiontab[-1,]
  reactant <- rownames(Mtotal)
  x_cond=as.matrix(x_cond)
  if (ncol(x_cond) == 1) {
    x_cond <- t(x_cond)
  }
  # Definitions ##################################################################
  
  # index for external reactants
 # n <- 1:dim(x_cond)[1]
  
  # internal matrix M
#  M <- Mtotal[-n,]
  
  # number of external reactants
 # nx <- dim(x_cond)[1]
  
  # number of growth conditions
  #n_conditions <- dim(x_cond)[2]
  
  # names of internal reactants
  #i_reactant <- reactant[-n]
  
  # number of internal reactants
  #p <- dim(M)[1]
  
  # number of reactions
  #r <- dim(M)[2]
  
  # the sum of each M column 
  #sM <- colSums(M)
  
  # delete numerical artifacts
  #sM[abs(sM) < 1e-10] <- 0
  
  # indexes for reactions: s (transport), e (enzymatic), and ribosome r 
  
  #e <- c(1:(r-1))[sM[1:(r-1)] == 0]  
  
  #s <- c(1:(r-1))[sM[1:(r-1)] != 0] 
  
  # indexes: m (metabolite), a (all proteins)
  
  #m <- 1:(p-1)
  
  # number of transport reactions
  #ns <- length(s)
  
  
  return(list(Mtotal = Mtotal, K = K, KI = KI, KA = KA, kcatf = kcatf, kcatb = kcatb,
              condition = condition, rho_cond = rho_cond, x_cond = x_cond,kcat_matrix=kcat_matrix,conditiontab=conditiontab))
}
