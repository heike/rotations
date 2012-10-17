
#' Create Euler angles
#' 
#' Create Euler angles representing the rotation of the identity matrix about the 
#' axis U throught the angle theta
#' 
#' @export
#' @param U three-dimensional vector describing the fix axis of the rotation
#' @param theta angle between -pi and pi
#' @return Euler angle of class "EA"

EA <- function(U,theta){
  #See if this can be made more efficient, rotations book
  R<-SO3(U,theta)
  eur<-EA.SO3(R)
  return(eur)
  
}


#'Translate from matrix to Euler angle format
#'
#'Based on the Z-X-Z definition of Euler angles, this will take a 3x3 rotation matrix and return
#'the Euler angle reprentation of that rotation.  See \cite{morawiec04} for a full description of this process.
#'@param rot a rotation matrix in the form of a 3x3 matrix
#'@return a vector of Euler angles
#'@cite morawiec04

EA.SO3 <- function(rot){
  
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

#' Identity in EA space
#' @export
id.EA <- as.EA(c(0,0,0))

#' Convert anything into EA class
#' 
#' @param x can be anything
#' @return x with class "EA"
#' @export

as.EA<-function(x){
  class(x)<-"EA"
  return(x)
}


#' Form a unit quaterion
#' 
#' Create a unit quaternion representing the rotation of the identity matrix about the 
#' axis U throught the angle theta
#' 
#' @export
#' @param U three-dimensional vector describing the fix axis of the rotation
#' @param theta angle between -pi and pi
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

Q4.SO3 <- function(eur) {
  
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
#' @param theta angle between -pi and pi
#' @return matrix of rotations in SO3 format

SO3 <- function(U, theta) {
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


#' A function that will take  a Euler angle and return a rotation matrix in vector format
#'
#' @param eur numeric Euler angle representation of an element in SO(3)
#' @return numeric \eqn{9\times 1} vector of a matrix in SO(3)
#' @seealso \code{\link{is.SO3}} can be used to check the output of this function
#' @export
#' @examples
#' eaExample<-structure(c(pi/2,3*pi/4,0), class="EA")
#' SO3Dat<-SO3.EA(eaExample)
#' is.SO3(SO3Dat)

SO3.EA <- function(eur) {
  
  e1<-c(1,0,0)
  e3<-c(0,0,1)
  S<-SO3(e3,eur[3])%*%SO3(e1,eur[2])%*%SO3(e3,eur[1])
  class(S)<-"SO3"
  return(S)
}

#' Translate a unit quaternion to a rotation matrix
#'
#' A function to translate from unit quaternion representation to \eqn{SO(3)} representation
#' of a rotation matrix.  Wikipedia has a good summary of this and other transforms.
#'
#' @param q numeric unit vector, i.e. \eqn{q^\top q=1}, representing an element in SO(3)
#' @return vector representation of a rotation matrix in SO(3)
#' @seealso \code{\link{is.SO3}} can be used to check the return vector
#' @export
#' @examples
#' is.SO3(SO3.Q4(c(1/sqrt(2),0,0,1/sqrt(2))))

SO3.Q4<-function(q){
  
  if((sum(q^2)-1)>10e-10){
    warning("Unit quaternions required.  Input was normalized.")
    q<-q/sqrt(sum(q^2))
  }
  
  
  theta<-angle(q)
  
  u<-axis2(q)
  
  return(SO3(u, theta)) 
}


#' Identity in SO(3) space
#' @export
id.SO3 <- as.SO3(diag(c(1,1,1)))

#' Convert anything into SO3 class
#' 
#' @param x matrix of rotations, note that no check is performed 
#' @return x with class "SO3"
#' @export

as.SO3<-function(x){
  class(x)<-"SO3"
  return(x)
}


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
  
  ## the second check will take care of the determinant  
  #  # Does it have determinant 1?
  #  if (abs(det(x)- 1)>10e-10) {
  #    return(FALSE)
  #  }
  
  # Is its transpose (approximately) its inverse?
  return(all(sum(t(x) %*% x - diag(1, 3))<10e-10)) 
  
}

