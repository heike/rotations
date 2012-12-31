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

#' Euler Angles
#' 
#' Create a Euler angles
#' 
#' Create Euler angles representing the rotation of the identity matrix about the 
#' axis U throught the angle theta.  This can be accomplished by providing the axis and angle
#' explicitly or by providing the rotation in some other form, e.g. a matrix in SO(3) or unit quaternion.
#' 
#' @export
#' @param U three-dimensional vector describing the fix axis of the rotation
#' @param ... additional arguments 
#' @return Euler angle of class "EA"
#' @family EA

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
#' @family EA

EA.SO3 <- function(R){  
  
	R<-formatSO3(R)
	trouble<-FALSE
  zeta<-sqrt(1-R[3,3]^2)
  
  #For now deal with singularity this way
  if(zeta==0){
    ea<-c(0,0,0)
    class(ea)<-"EA"
    return(ea)
  }
  
  Salpha <- asin(R[3,1]/zeta)
  if(Salpha<0)  Salpha<-Salpha+2*pi
  
  Calpha <- acos(-R[3,2]/zeta)
  if(Calpha<0)  Calpha<-Calpha+2*pi
  
  if(Salpha==Calpha){
    alpha <- Salpha 
  }else if((sin(Calpha)-(R[3,1]/zeta))<10e-5){
    alpha <- Calpha
  }else if((cos(Salpha)+R[3,2]/zeta)<10e-5){
    alpha <- Salpha
  }else{
    alpha<-Salpha
    trouble<-TRUE
  }
  
  beta<-acos(R[3,3])
  
  Sgamma <- asin(R[1,3]/zeta)
  if(Sgamma<0) Sgamma<-Sgamma+2*pi
  
  Cgamma<- acos(R[2,3]/zeta)
  if(Cgamma<0) Cgamma<-Cgamma+2*pi
  
  if(Sgamma==Cgamma){
    gamma <- Sgamma 
  }else if((sin(Cgamma)-R[1,3]/zeta)<10e-5){
    gamma <- Cgamma    
  }else if((cos(Sgamma)-R[2,3]/zeta)<10e-5){
    gamma <- Sgamma
  }else{
    gamma<-Sgamma
    trouble<-TRUE
  }
    
  ea<-c(alpha,beta,gamma)
	
	if(trouble)
    ea[c(1,3)]<-(ea[c(1,3)]+pi)%%(2*pi)
  
	#if(any(ea<0)){
	#	ea[c(1,3)]<-ea[c(1,3)]+pi
	#	ea[2]<-2*pi-ea[2]
	#}
	
  class(ea)<-"EA"
  return(ea)
}

#' @return \code{NULL}
#' 
#' @rdname EA
#' @method EA Q4
#' @S3method EA Q4
#' @family EA

EA.Q4 <- function(Qs){
  theta<-angle(Qs)
  u<-axis2(Qs)
  return(EA(u,theta))
}

#' Convert anything into EA class
#' 
#' @param x can be anything
#' @return x with class "EA"
#' @family EA
#' @export

as.EA<-function(x){
  class(x)<-"EA"
  return(x)
}

#' Identity in EA space
#' @family EA
#' @export
id.EA <- as.EA(c(0,0,0))


#' Quaternions
#' 
#' Create a unit quaternion
#' 
#' Create quaternion representing the rotation of the identity matrix about the 
#' axis U throught the angle theta.  This can be accomplished by providing the axis and angle
#' explicitly or by providing the rotation in some other form, e.g. a matrix in SO(3) or Euler angles.
#' 
#' @export
#' @param U three-dimensional vector describing the fix axis of the rotation
#' @param ... additional arguments
#' @return unit quaternion of class "Q4"
#' @family Q4

Q4<-function(U,...){
  UseMethod("Q4")
}

#' @return \code{NULL}
#'
#' @rdname Q4
#' @method Q4 default
#' @S3method Q4 default
#' @family Q4

Q4.default <- function(U,theta){
	n<-length(theta)
	
	if(n!=length(U)/3){
		stop("A different number of angles and axis were supplied.")
	}
	
	x <- cbind(cos(theta/2), sin(theta/2) * U)

	class(x)<-"Q4"
  return(x)
}

#' @return \code{NULL}
#'
#' @rdname Q4
#' @method Q4 SO3
#' @S3method Q4 SO3
#' @family Q4

Q4.SO3 <- function(R) {
  
	R<-formatSO3(R)
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
#' @family Q4

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
#' @family Q4
#' @export

as.Q4<-function(x){
  class(x)<-"Q4"
  return(x)
}

#' Identity in Q4 space
#' @family Q4
#' @export
id.Q4 <- as.Q4(c(1,0,0,0))

#' Matrix in SO(3)
#' 
#' Create a rotation matrix
#' 
#' Create matrix in SO(3) representing the rotation of the identity matrix about the 
#' axis U throught the angle theta.  This can be accomplished by providing the axis and angle
#' explicitly or by providing the rotation in some other form, e.g. a vector of Euler angles or unit quaternion.
#'
#' @export
#' @param U three-dimensional vector describing the fix axis of the rotation
#' @param ... additional arguments
#' @return matrix of rotations in SO3 format
#' @family SO3

SO3 <- function(U,...){
  UseMethod("SO3")
}

#' @return \code{NULL}
#' 
#' @rdname SO3
#' @method SO3 default
#' @S3method SO3 default
#' @family SO3

SO3.default <- function(U, theta=NULL) {
  
	n<-length(U)/3
	
	if(n%%1!=0)
		stop("This functions only works in three dimensions.")	
	
	U<-matrix(U,n,3)
	
	ulen<-sqrt(rowSums(U^2)) 
  
  if(is.null(theta)){ 
  	theta<-ulen%%(2*pi)
  	
  	#if(theta>pi)
  	#	theta<-2*pi-theta
  }

	R<-matrix(NA,n,9)
	
  for(i in 1:n){
  	
  	if(ulen[i]!=0)
  		U[i,]<-U[i,]/ulen[i]
  	
 		P <- U[i,] %*% t(U[i,])
   
  	R[i,] <- P + (diag(3) - P) * cos(theta[i]) + eskew(U[i,]) * sin(theta[i])
  }
 		
 	class(R) <- "SO3"
  return(R)
}


#' @return \code{NULL}
#' 
#' @rdname SO3
#' @method SO3 EA
#' @S3method SO3 EA
#' @family SO3

SO3.EA <- function(eur) {
  
  S1<-S2<-S3<-diag(3)
  S1[1,1]<-cos(eur[3])
  S1[1,2]<-sin(eur[3])
  S1[2,1]<--sin(eur[3])
  S1[2,2]<-cos(eur[3])
  
  S2[2,2]<-cos(eur[2])
  S2[2,3]<-sin(eur[2])
  S2[3,2]<--sin(eur[2])
  S2[3,3]<-cos(eur[2])
  
  S3[1,1]<-cos(eur[1])
  S3[1,2]<-sin(eur[1])
  S3[2,1]<--sin(eur[1])
  S3[2,2]<-cos(eur[1])
  
  S<-S1%*%S2%*%S3
  class(S)<-"SO3"
  return(S)
}

#' @return \code{NULL}
#' 
#' @rdname SO3
#' @method SO3 Q4
#' @S3method SO3 Q4
#' @family SO3

SO3.Q4<-function(q){
  
  if(any((rowSums(q^2)-1)>10e-10)){
    warning("Unit quaternions required.  Input was normalized.")
    nonq<-which((rowSums(q^2)-1)>10e-10)
    q[nonq,]<-as.Q4(q[nonq,]/sqrt(rowSums(q[nonq,]^2)))
  }else{
  	q<-as.Q4(q)
  }
  
  theta<-angle(q)
  
  u<-axis2(q)
  
  return(SO3(u, theta)) 
}


#' Convert anything into SO3 class
#' 
#' @param x matrix of rotations, note that no check is performed 
#' @return x with class "SO3"
#' @family SO3
#' @export

as.SO3<-function(x){
  class(x)<-"SO3"
  return(x)
}

#' Identity in SO(3) space
#' @family SO3
#' @export
id.SO3 <- as.SO3(diag(c(1,1,1)))


#' A function to determine if a given matrix is in \eqn{SO(3)} or not.
#'
#' @param x numeric 3-by-3 matrix or vector of length 9
#' @return logical T if the matrix is in SO(3) and false otherwise
#' @family SO3
#' @export
#' @examples
#' is.SO3(diag(1,3,3))
#' is.SO3(1:9)
is.SO3 <- function(x) {
  
  x <- matrix(x, 3, 3)
  if (any(is.na(x))) return(FALSE)
  
  return(all(sum(t(x) %*% x - diag(1, 3))<10e-10)) 
  
}

