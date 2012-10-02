#'  Find kappa for given nu
#'  
#'  @param nu The circular variance
#'  @return the concentration parameter corresponding to nu
#'  @export

cayley_kappa<-function(nu){
  (3/nu)-2
}

fisher_nu_kappa<-function(kappa,nu){
  (1-(besselI(2*kappa,1)-.5*besselI(2*kappa,2)-.5*besselI(2*kappa,0))/(besselI(2*kappa,0)-besselI(2*kappa,1))-nu)^2
}

#'  Find kappa for given nu
#'  
#'  @param nu The circular variance
#'  @return the concentration parameter corresponding to nu
#'  @export

fisher_kappa<-function(nu){
  
  kappa<-rep(0,length(nu))
  
  for(i in 1:length(nu))
    kappa[i]<-optimize(fisher_nu_kappa,interval=c(0,10),tol=.00001,nu=nu[i])$minimum
  
  return(kappa)
}


mises_nu_kappa<-function(kappa,nu){
  (1-besselI(kappa,1)/besselI(kappa,0)-nu)^2
}

#'  Find kappa for given nu
#'  
#'  @param nu The circular variance
#'  @return the concentration parameter corresponding to nu
#'  @export

vmises_kappa<-function(nu){
  
  kappa<-rep(0,length(nu))
  
  for(i in 1:length(nu))
    kappa[i]<-optimize(mises_nu_kappa,interval=c(0,10),tol=.00001,nu=nu[i])$minimum
  
  return(kappa)
}
