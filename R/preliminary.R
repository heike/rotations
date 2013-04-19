#' Accept/Reject Algorithm
#'
#' @author Heike Hofmann
#' @param f target density
#' @param g sampling density
#' @param M real valued constant, maximum of g
#' @param kappa second parameter in the target density
#' @param Haar should f be with respect to Haar or not
#' @param ... additional arguments passed to samping density, g
#' @return a random observation from target density

arsample <- function(f, g, M, kappa, Haar, ...) {
  found = FALSE
  while (!found) {
    x <- g(1, ...)
    y <- runif(1, min = 0, max = M)
    if (y < f(x, kappa, Haar)) 
      found = TRUE
  }
  return(x)
  # arsample(f, g, M, kappa, ...)
}

#' Accept/Reject Algorithm
#'
#' @author Heike Hofmann
#' @param f target density
#' @param M maximum value for enveloping uniform density
#' @param ... additional arguments sent to f
#' @return x an observation from the target density

arsample.unif <- function(f, M, ...) {
  found = FALSE
  while (!found) {
    x <- runif(1, -pi, pi)
    y <- runif(1, min = 0, max = M)
    if (y < f(x, ...)) 
      found = TRUE
  }
  return(x)
  # arsample.unif(f, M, ...)
}



#' Rotational Distance
#'
#' Calculate the Euclidean or Riemannian distance between two rotations
#'
#' This function will calculate the intrinsic (Riemannian) or projected (Euclidean) distance between two rotations.  If only one rotation is specified
#' the other will be set to the identity and the distance between the two is returned.  For rotations \eqn{R_1}{R1} and \eqn{R_2}{R2}
#' both in \eqn{SO(3)}, the Euclidean distance between them is \deqn{||R_1-R_2||_F}{||R1-R2||} where \eqn{||\cdot||_F}{|| ||} is the Frobenius norm.
#' The intrinsic distance is defined as \deqn{||Log(R_1^\top R_2)||_F}{||Log(R1'R2)||} where \eqn{Log} is the matrix logarithm, and it corresponds
#' to the misorientation angle of \eqn{R_1^\top R_2}{R1'R2}.
#'
#' @param R1 (or Q1) a rotation in matrix or quaternion representation
#' @param R2 (or Q2) the second rotation in the same parameterization as R1
#' @param method String indicating 'projected' or 'intrinsic' method of distance 
#' @param p the order of the distance 
#' @param ... Additional arguments
#' @return the rotational distance between R1 and R2
#' @export

dist<-function(x,...){
  UseMethod("dist")
}


#' @rdname dist
#' @method dist SO3
#' @S3method dist SO3

dist.SO3 <- function(R1, R2=id.SO3, method='projected' , p=1) {
  
  R1<-formatSO3(R1)
  
  if(method=='projected'){
    
    R2<-matrix(R2,nrow(R1),9,byrow=T)
    so3dist<-sqrt(rowSums((R1-R2)^2))^p
    
  }else if(method=='intrinsic'){
    
    R2<-matrix(R2,3,3)
    
    for(i in 1:nrow(R1)){
      R1[i,]<-t(matrix(R1[i,],3,3))%*%R2
    }
    
    so3dist<-angle(as.SO3(R1))^p
    
  }else{
    stop("Incorrect usage of method argument.  Please choose intrinsic or projected.")
  }
  
  return(so3dist)
  
}


#' @rdname dist
#' @method dist Q4
#' @S3method dist Q4

dist.Q4 <- function(Q1, Q2=id.Q4 ,method='projected', p=1) {

  Q1<-formatQ4(Q1)
  Q2<-formatQ4(Q2)
  
  if(method=='intrinsic'){
    
    Q2<-matrix(Q2,nrow(Q1),4,byrow=T)
    cp <- rowSums(Q1*Q2)
    q4dist<-acos(2*cp*cp-1)^p
    
  }else if(method=='projected'){
    
    R1<-SO3.Q4(Q1)
    R2<-SO3.Q4(Q2)
    R2<-matrix(R2,nrow(R1),9,byrow=T)
    q4dist<-sqrt(rowSums((R1-R2)^2))^p
    
  }else{
    stop("Incorrect usage of method argument.  Please choose intrinsic or projected.")
  }
  
  return(q4dist)
}


#' Misorientation Angle
#' 
#' Find the misorientation angle of a rotation
#' 
#' Every rotation can be thought of as some reference coordinate system rotated about an axis through an angle.  These quantites
#' are referred to as the misorientation axis and misorientation angle, respectively, in the material sciences literature.
#' This function returns the misorentation angle associated with a rotation assuming the reference coordinate system
#' is the identity.
#'  
#' @param Rs rotation matrix
#' @return angle of rotation
#' @seealso \code{\link{axis}}
#' @export

angle<-function(Rs){
  UseMethod("angle")
}


#' @rdname angle
#' @method angle SO3
#' @S3method angle SO3

angle.SO3 <- function(Rs){
	
	Rs<-formatSO3(Rs)
	n<-nrow(Rs)
	tr<-rep(0,n)
	eps <- 10^-3
	for(i in 1:n){
  	##  trace of a rotation matrix has to be between -1 and 3. If not, this is due
  	## to numerical inconcistencies, that we have to fix here
  	tr[i]<-Rs[i,1]+Rs[i,5]+Rs[i,9]

		if (tr[i] > 3+eps) print(sprintf("Warning: trace too large (> 3.0): %f ", tr))
  		if (tr[i] < -1-eps) print(sprintf("Warning: trace too small (< -1.0): %f ", tr))
		#  stopifnot(tr<3+eps, tr>-1-eps)
  	tr[i] <- max(min(3, tr[i]), -1)
	}
  return(acos((tr-1)/2))
}


#' @rdname angle
#' @method angle Q4
#' @S3method angle Q4

angle.Q4 <- function(Qs){
	
  Qs<-formatQ4(Qs)
	n<-nrow(Qs)
	theta<-rep(0,n)
	
	for(i in 1:n){
	  theta[i]<-2*acos(Qs[i,1])
	}
	 return(theta)
}


#' Misorientation Axis
#' 
#' Find the misorientation axis of a rotation
#' 
#' Every rotation can be thought of as some reference coordinate system rotated about an axis through an angle.  These quantites
#' are referred to as the misorientation axis and misorientation angle, respectively, in the material sciences literature.
#' This function returns the misorentation axis associated with a rotation assuming the reference coordinate system
#' is the identity.
#' 
#' @param R 3-by-3 matrix in SO3 
#' @return axis in form of three dimensional vector of length one.
#' @seealso \code{\link{angle}}
#' @export

axis2<-function(R){
  UseMethod("axis2")
}

#' @rdname axis2
#' @method axis2 SO3
#' @S3method axis2 SO3

axis2.SO3<-function(R){
  
	R<-formatSO3(R)
  n<-nrow(R)
	u<-matrix(NA,n,3)
	
	for(i in 1:n){
		Ri<-matrix(R[i,],3,3)
  	X <- Ri - t(Ri)
  	u[i,] <- rev(X[upper.tri(X)])*c(-1,1,-1)
		u[i,]<-u[i,]/sqrt(sum(u[i,]^2))
	}
  return(u) # will be trouble, if R is symmetric, i.e. id,  .... 

}

#' @rdname axis2
#' @method axis2 Q4
#' @S3method axis2 Q4

axis2.Q4 <- function(q){
  
  q<-formatQ4(q)
  theta<-angle(q)
  
  u <- q[,2:4]/sin(theta/2)

	if(any(is.infinite(u)|is.nan(u))){
		infs<-which(is.infinite(u)|is.nan(u))
		u[infs]<-0
	}  
  
  return(u)
}


#' Directional vector to skew-symmetric Matrix
#'
#' @author Heike Hofmann
#' @param U three dimensional vector indicating rotational fix-axis
#' @return skew-symmetric matrix

eskew <- function(U) {
  
  ulen<-sqrt(sum(U^2))
  
  if(ulen!=0){
    U<-U/ulen
  }
  
  u <- U[1]
  v <- U[2]
  w <- U[3]
  
  res <- matrix((-1) * c(0, -w, v, w, 0, -u, -v, u, 0), ncol = 3)
  return(res)
}



#' Generate Rotations
#'
#' Generate rotations according to the Uniform-Axis Random Spin methodology
#'
#' Given a vector \eqn{u\in\mathbb{R}^2}{u in R^2} of length one and angle of rotation \eqn{r}, a rotation can be formed using Rodrigues formula
#' \deqn{\cos(r)I_{3\times 3}+\sin(r)\Phi(u)+(1-\cos(r))uu^\top}{cos(r)I+sin(r)\Phi(u)+(1-cos(r))uu'} 
#' where \eqn{I_{3\times 3}}{I} is the \eqn{3\times 3}{3-by-3} identity matrix,\eqn{\Phi(u)} is a \eqn{3\times 3}{3-by-3} skew-symmetric matirix
#' with upper triangular elements \eqn{-u_3}{-u3}, \eqn{u_2}{u2} and \eqn{-u_1}{-u1} in that order.
#'
#' @param r vector of angles
#' @param S The principle direction
#' @param space Indicates the desired representation: matrix in SO3, quaternion, or Euler angles 
#' @return a matrix where each row is a sample point in the desired space
#' @cite bingham09
#' @export
#' @examples
#' r<-rvmises(20,0.01)
#' genR(r)

genR <- function(r, S = NULL, space='SO3') {
  
  if(!(space %in% c("SO3","Q4")))
    stop("Incorrect space argument.  Options are: SO3 and Q4. ")
  
  n<-length(r)
  
  # Generate angles theta from a uniform distribution from 0 to pi
  
  theta <- acos(runif(n, -1, 1))
  
  # Generate angles phi from a uniform distribution from -pi to pi
  
  phi <- runif(n, -pi, pi)
  u <- matrix(c(sin(theta) * cos(phi), sin(theta) * sin(phi), cos(theta)),length(r),3)
  
  if(space=="SO3"){
  	
  	if(is.null(S))
  		S<-id.SO3
  	S<-formatSO3(S)
  	o<-SO3(u,r)
  	o<-centeringSO3(o,t(S))
  	
  	class(o) <- "SO3"
  	return(o)
  	
  }else{
  	
  	if(is.null(S))
  		S<-id.Q4
  	
  	S<-formatQ4(S)
  	
  	S[2:4]<--S[2:4]
  	q<-matrix(c(cos(r/2),sin(r/2)*u),n,4)
  	q<-centeringQ4(q,S)
  	
  	class(q)<-"Q4"
  	return(q)
  }

}


#' This fuction will compute the natural exponential of skew-symmetric matrix.
#'
#' See \cite{moakher02}
#'
#' @param A 3-dimensional skew-symmetric matrix, i.e., \eqn{\bm A=-\bm A^\top}
#' @return numeric matrix \eqn{e^{\bm A}}{e^A}
#' @cite moakher02

exp.skew <- function(A) {
  
  if (sum(abs(A + t(A)))>10e-10) {
    stop("The input matrix must be skew symmetric.")
  }
  
  I <- diag(1, 3, 3)
  AtA <- t(A) %*% A
  a2 <- 0.5 * sum(diag(AtA))
  a <- sqrt(a2)
  
  if (a == 0) {
    return(I)
  } else {
    p1 <- (sin(a)/a) * A
    p2 <- ((1 - cos(a))/a^2) * A %*% A
    return(I + p1 + p2)
  }
  
}


#' Natural Logarithm in SO(3)
#'
#' For details see \cite{moakher02}
#'
#' @param R numeric matrix in \eqn{SO(n)}
#' @return mlog numeric matrix \eqn{\log(R)}{log(R)}
#' @cite moakher02

log.SO3 <- function(R) {
  
	#R<-formatSO3(R)
  
  theta <- angle.SO3(R)  
  R<-matrix(R,3,3)
  
  if (abs(cos(theta)) >= 1) {
    return(diag(0, 3, 3))
  } else {
    
    c <- theta/(2 * sin(theta))
    mlog <- c * (R - t(R))
    return(mlog)
  }
}

#' Projection Procedure
#'
#' Project an arbitrary \eqn{3\times 3}{3-by-3} matrix into SO(3)
#'
#' This function uses the process given in \cite{moakher02} to project an arbitrary \eqn{3\times 3}{3-by-3} matrix into \eqn{SO(3)}.
#' 
#' @param M \eqn{3\times 3}{3-by-3} matrix to project
#' @return projection of \eqn{\bm M}{M} into \eqn{SO(3)}
#' @seealso \code{\link{mean.SO3}}, \code{\link{median.SO3}}
#' @export
#' @examples
#' M<-matrix(rnorm(9),3,3)
#' project.SO3(M)

project.SO3 <- function(M) {
  
	M<-matrix(M,3,3)
  d <- svd(t(M) %*% M)
  u <- d$u
  
  d1 <- 1/sqrt(d$d[1])
  d2 <- 1/sqrt(d$d[2])
  d3 <- sign(det(M))/sqrt(d$d[3])
  
  R <- M %*% u %*% diag(x = c(d1, d2, d3), 3, 3) %*% t(u)
  return(R)
}


#' Sample Distance
#'
#' Compute the sum of the \eqn{p^{\text{th}}}{pth} order distances between Rs and S
#'
#' @param Rs a matrix of rotation observations, one row per observation
#' @param S the individual matrix of interest, usually an estimate of the mean
#' @param method type of distance used method in 'projected' or 'intrinsic'
#' @param p the order of the distances to compute
#' @return the sum of the pth order distance between each sample in Rs and S
#' @export
#' @examples
#' r<-rvmises(20,0.01)
#' Rs<-genR(r)
#' Sp<-mean(Rs)
#' sum_dist(Rs,S=Sp,p=2)

sum_dist<-function(Rs, S = genR(0, space=class(Rs)), method='projected', p=1){
  
  UseMethod( "sum_dist" )

}

#' @return \code{NULL}
#'
#' @rdname sum_dist
#' @method sum_dist SO3
#' @S3method sum_dist SO3

sum_dist.SO3 <- function(Rs, S = id.SO3, method='projected', p=1) {

  return(sum(dist(Rs,S, method=method, p=p)))
  
}


#' @return \code{NULL}
#'
#' @rdname sum_dist
#' @method sum_dist Q4
#' @S3method sum_dist Q4

sum_dist.Q4 <- function(Qs, S = id.Q4, method='projected', p=1) {
  
  return(sum(dist(Qs,S, method=method, p=p)))
  
}



tLogMat <- function(x, S) {
  tra <- log.SO3(t(S) %*% matrix(x, 3, 3))
  return(as.vector(tra))
}

tLogMat2 <- function(x, S) {
  tra <- log.SO3(matrix(x, 3, 3)%*%t(S))
  return(as.vector(tra))
}

vecNorm <- function(x, S, ...) {
  n <- sqrt(length(x))
  cenX <- x - as.vector(S)
  return(norm(matrix(cenX, n, n), ...))
}


#' Centering function
#' 
#' This function will take the sample Rs and return teh sample Rs centered at
#' S, i.e., the returned sample is S'Rs, so if S is the true center then
#' the projected mean should be id.SO3
#' 
#' @param Rs the sample to be centered
#' @param S the rotation to sample around
#' @return The centered sample
#' @export

centeringSO3<-function(Rs,S){
	#This takes a set of observations in SO3 and centers them around S
	
	Rs<-formatSO3(Rs)
	S<-matrix(formatSO3(S),3,3)
	
	for(i in 1:nrow(Rs)){
		Rs[i,]<-t(S)%*%matrix(Rs[i,],3,3)
	}
	return(as.SO3(Rs))
}

centeringQ4<-function(Qs,S){
	#This takes a set of observations in Q4 and centers them around S
	Qs<-formatQ4(Qs)
	S<-formatQ4(S)
	S[2:4]<--S[2:4]
	
	for(i in 1:nrow(Qs)){
		Qs[i,]<-qMult(S,Qs[i,])
	}

	return(Qs)
}

formatSO3<-function(Rs){
	#This function will take input and format it to work with our functions
	#It also checks that the data is actually SO3 and of appropriate dimension
	
	len<-length(Rs)
	if(len%%9!=0)
		stop("Data needs to have length divisible by 9.")
	
	Rs<-matrix(Rs,len/9,9)
	
	if (!all(apply(Rs, 1, is.SO3))) 
		warning("At least one of the given observations is not in SO(3).  Use result with caution.")
	
	
	return(as.SO3(Rs))

}

formatQ4<-function(Qs){
	
  if(length(Qs)%%4!=0)
    stop("Data needs to have length divisible by 4.")
  
  Qs<-matrix(Qs,length(Qs)/4,4)
  
  if (!all(apply(Qs, 1, is.Q4))) 
  	warning("At least one of the given observations is not a unit quaternion.  Use result with caution.")
  
  
  if(length(Qs)==4)
    return(as.Q4(Qs))
  else
    return(as.Q4(Qs))
}

pMat<-function(p){
	#Make the matrix P from quaternion p according to 3.1 of Rancourt, Rivest and Asselin (2000)
	#This is one way to multiply quaternions
	Pmat<-matrix(0,4,4)
	Pmat[,1]<-p
	Pmat[,2]<-p[c(2,1,4,3)]*c(-1,1,1,-1)
	Pmat[,3]<-c(-p[3:4],p[1:2])
	Pmat[,4]<-p[4:1]*c(-1,1,-1,1)
	return(Pmat)
}

qMult<-function(q1,q2){
	#Forms quaternion product q1 x q2, i.e., rotate q2 by q1
	#This functions utilizes the 
	q1<-formatQ4(q1)
	q2<-formatQ4(q2)
	q1q2<-pMat(q1)%*%matrix(q2,4,1)
	return(formatQ4(q1q2))
}

proj<-function(u,v){
	#Project the vector v orthogonally onto the line spanned by the vector u
	num<-t(u)%*%v
	denom<-t(u)%*%u
	return(num*u/denom)
}
