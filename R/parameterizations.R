#' "SO3" class
#'
#' @name SO3-class
#' @aliases SO3
#' @family SO3
#'
#' @exportClass SO3
setOldClass("SO3")


#' Q4 class
#'
#' Class for quaterion representation of rotations
#' 
#' @name Q4-class
#' @aliases Q4
#' @family Q4
#'
#' @exportClass Q4
setOldClass("Q4")

#' EA class
#'
#' Class for Euler angle representation of rotations
#' 
#' @name EA-class
#' @aliases EA
#' @family EA
#'
#' @exportClass EA
setOldClass("EA")

#' Create Euler angles
#' 
#' Create Euler angles representing the rotation of the identity matrix about the 
#' axis U throught the angle theta
#' 
#' @export
#' @param U three-dimensional vector describing the fix axis of the rotation
#' @param ... additional arguments 
#' @return Euler angle of class "EA"

EA<-function(U,...){
  UseMethod("EA")
}

#' @return \code{NULL}
#' 
#' @rdname EA
#' @method EA default
#' @S3method EA default

EA.default <- function(U,theta){
  #See if this can be made more efficient, rotations book
  R<-SO3(U,theta)
  eur<-EA.SO3(R)
  return(eur)
  
}


#' @return \code{NULL}
#' 
#' @rdname EA
#' @method EA SO3
#' @S3method EA SO3

EA.SO3 <- function(R){  
  
  zeta<-sqrt(1-rot[9]^2)
  
  #For now deal with singularity this way
  if(zeta==0){
    ea<-c(0,0,0)
    class(ea)<-"EA"
    return(ea)
  }
  
  Salpha <- asin(rot[3]/zeta)
  if(Salpha<0)  Salpha<-Salpha+2*pi
  
  Calpha <- acos(-rot[6]/zeta)
  
  if(Calpha<0)  Calpha<-Calpha+2*pi
  
  if(Salpha==Calpha){
    alpha <- Salpha 
  }else if(sin(Calpha)==rot[3]/zeta){
    alpha <- Calpha
  }else{
    alpha <- Salpha
  }
  
  
  beta<-acos(rot[9])
  
  Sgamma <- asin(rot[7]/zeta)
  if(Sgamma<0) Sgamma<-Sgamma+2*pi
  
  Cgamma<- acos(rot[8]/zeta)
  if(Cgamma<0) Cgamma<-Cgamma+2*pi
  
  if(Sgamma==Cgamma){
    gamma <- Sgamma 
  }else if(sin(Cgamma)==rot[3]/zeta){
    gamma <- Cgamma    
  }else{
    gamma <- Sgamma
    
  }
  
  ea<-c(alpha,beta,gamma)
  class(ea)<-"EA"
  return(ea)
}

#' @return \code{NULL}
#' 
#' @rdname EA
#' @method EA Q4
#' @S3method EA Q4

EA.Q4 <- function(Qs){
  theta<-angle(Qs)
  u<-axis2(Qs)
  return(EA(u,theta))
}

#' Convert anything into EA class
#' 
#' @param x can be anything
#' @return x with class "EA"
#' @export

as.EA<-function(x){
  class(x)<-"EA"
  return(x)
}

#' Identity in EA space
#' @export
id.EA <- as.EA(c(0,0,0))


#' Form a unit quaterion
#' 
#' Create a unit quaternion representing the rotation of the identity matrix about the 
#' axis U throught the angle theta
#' 
#' @export
#' @param U three-dimensional vector describing the fix axis of the rotation
#' @param ... additional arguments
#' @return unit quaternion of class "Q4"

Q4<-function(U,...){
  UseMethod("Q4")
}

#' @return \code{NULL}
#'
#' @rdname Q4
#' @method Q4 default
#' @S3method Q4 default

Q4.default <- function(U,theta){
  x <- c(cos(theta/2), sin(theta/2) * U)
  class(x)<-"Q4"
  return(x)
}

#' @return \code{NULL}
#'
#' @rdname Q4
#' @method Q4 SO3
#' @S3method Q4 SO3

Q4.SO3 <- function(R) {
  
  theta <- angle(R)
  u <- axis2(R)
  x <- Q4(u,theta)
  
  return(x)
}

#' @return \code{NULL}
#'
#' @rdname Q4
#' @method Q4 EA
#' @S3method Q4 EA

Q4.EA <- function(eur) {
  
  theta <- angle(eur)
  u <- axis2(eur)
  x <- Q4(u,theta)
  
  return(x)
}

#' Convert anything into Q4 class
#' 
#' @param x can be anything
#' @return x with class "Q4"
#' @export

as.Q4<-function(x){
  class(x)<-"Q4"
  return(x)
}

#' Identity in Q4 space
#' @export
id.Q4 <- as.Q4(c(1,0,0,0))

#' Method for creating a rotation in SO3 format using the angle axis representation
#'
#' Angle-axis representation based on the Rodrigues formula.
#'
#' @export
#' @param U three-dimensional vector describing the fix axis of the rotation
#' @param ... additional arguments
#' @return matrix of rotations in SO3 format

SO3 <- function(U,...){
  UseMethod("SO3")
}

#' @return \code{NULL}
#' 
#' @rdname SO3
#' @method SO3 default
#' @S3method SO3 default

SO3.default <- function(U, theta) {
  # based on Rodrigues formula
  
  ulen<-sqrt(sum(U^2))
  
  if(ulen!=0){
    U <- U/ulen
  }
  
  P <- U %*% t(U)
  
  id <- matrix(0, length(U), length(U))
  diag(id) <- 1
  
  
  R <- P + (id - P) * cos(theta) + eskew(U) * sin(theta)
  class(R) <- "SO3"
  return(R)
}


#' @return \code{NULL}
#' 
#' @rdname SO3
#' @method SO3 EA
#' @S3method SO3 EA

SO3.EA <- function(eur) {
  
  e1<-c(1,0,0)
  e3<-c(0,0,1)
  S<-SO3(e3,eur[3])%*%SO3(e1,eur[2])%*%SO3(e3,eur[1])
  class(S)<-"SO3"
  return(S)
}

#' @return \code{NULL}
#' 
#' @rdname SO3
#' @method SO3 Q4
#' @S3method SO3 Q4

SO3.Q4<-function(q){
  
  if((sum(q^2)-1)>10e-10){
    warning("Unit quaternions required.  Input was normalized.")
    q<-as.Q4(q/sqrt(sum(q^2)))
  }
  
  theta<-angle(q)
  
  u<-axis2(q)
  
  return(SO3(u, theta)) 
}


#' Convert anything into SO3 class
#' 
#' @param x matrix of rotations, note that no check is performed 
#' @return x with class "SO3"
#' @export

as.SO3<-function(x){
  class(x)<-"SO3"
  return(x)
}

#' Identity in SO(3) space
#' @export
id.SO3 <- as.SO3(diag(c(1,1,1)))


#' A function to determine if a given matrix is in \eqn{SO(3)} or not.
#'
#' @param x numeric \eqn{n \times n} matrix or vector of length \eqn{n^2}
#' @return logical T if the matrix is in SO(3) and false otherwise
#' @export
#' @examples
#' is.SO3(diag(1,3,3))
#' is.SO3(1:9)
is.SO3 <- function(x) {
  
  x <- matrix(x, 3, 3)
  if (any(is.na(x))) return(FALSE)
  
  return(all(sum(t(x) %*% x - diag(1, 3))<10e-10)) 
  
}
