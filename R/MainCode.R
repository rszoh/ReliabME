#----------------------------------------------------
#---- R code that accompany the paper titled: 
#---- Adjusting for Covariates Representing Potential Confounders, Mediators,
#---- or Competing Predictors in the Presence of Measurement Error: 
#---- Dispelling a Potential Misapprehension and Insights for Optimal Study Design with Nutritional Epidemiology Examples
#---- The paper was submitted to F1000
#----------------------------------------------------
require(numDeriv)
require(latex2exp)


#--- Partial correlation Y.X|Z' between Y an X, conditional on Z' (R^{2}_{Y,X.Z'})
ParCorsq1 <- function(gam,beta, alphaz, w){
  w = 1
  alphax = 1 
  
  Num <- (((alphax*gam*beta)^2 + ((alphaz*beta)^2)  - 2*((alphax*beta*alphaz*gam)^2) )/(1 - (alphax*gam*alphaz)^2)) - (alphaz*beta)^2
  Den  <- 1 - (alphaz*beta)^2
  
  Num/Den
}


#--- Partial correlation between Y an X', conditional on Z' (R^{2}_{Y,X'.Z'})
ParCorsq0 <- function(gam,beta, alphaz, w, a0){
  alphax = a0 + alphaz*w
  
  Num <- (((alphax*gam*beta)^2 + ((alphaz*beta)^2)  - 2*((alphax*beta*alphaz*gam)^2) )/(1 - (alphax*gam*alphaz)^2)) - (alphaz*beta)^2
  Den  <- 1 - (alphaz*beta)^2
  
  Num/Den
}
ParCorsq0 <- Vectorize(ParCorsq0,"alphaz") 

#--- This function will compute tghe partial derivative with respect to alphaz

DevFun <- function(alphaz, gam,beta,w){
  alphax = alphaz*w
  
  Num <- ((alphaz*beta*w*gam)^2)*(1 - 2*(alphaz^2) + (alphaz^4))
  Den = (1 - (alphaz*beta)^2)*(1 - (w*gam*(alphaz^2))^2)
  
  #Num <- (((alphax*gam*beta)^2 + (alphaz*beta)^2  - 2*(alphax*beta*alphaz*gam)^2)/(1 - (alphax*gam*alphaz)^2)) - (alphaz*beta)^2
  #Den  <- 1 - (alphaz*beta)^2
  Num/Den
  
}

Dev0 <- grad(DevFun, alphaz, gam = gam, beta=beta, w=-1)

