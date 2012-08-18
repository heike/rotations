#' Method for creating a rotation using the angle axis representation
#'
#' Angle-axis representation based on the Rodrigues formula.
#'
#' @export
#' @param U three-dimensional vector describing the fix axis of the rotation
#' @param theta angle between -pi and pi
#' @return Used in \code{\link{eyeBall}} to orient the data properly

angle_axis <- function(U, theta) {
  # based on Rodrigues formula
  U <- U/sqrt(sum(U^2))
  P <- U %*% t(U)
  
  id <- matrix(0, length(U), length(U))
  diag(id) <- 1
  
  
  R <- P + (id - P) * cos(theta) + eskew(U) * sin(theta)
  return(R)
}


#' Accept/reject algorithm random sampling from angle distributions
#' 
#' @author Heike Hofmann
#' @param f target density
#' @param g sampling density
#' @param M real valued constant, maximum of g
#' @param kappa second parameter in the target density
#' @param ... additional arguments passed to samping density, g
#' @return a random observation from target density

arsample <- function(f,g,M, kappa, ...) {
  found=FALSE
  while(!found) {
    x <- g(1, ...)
    y <- runif(1, min=0, max=M)
    if (y < f(x, kappa)) found=TRUE
  }
  return(x)
  #  arsample(f, g, M, kappa, ...)
}

#' Accept/reject algorithm random sampling from angle distributions using uniform envelop  
#' 
#' @author Heike Hofmann
#' @param f target density
#' @param M maximum value for enveloping uniform density
#' @param ... additional arguments sent to f
#' @return x an observation from the target density

arsample.unif <- function(f,M, ...) {
  found=FALSE
  while(!found) {
    x <- runif(1, -pi, pi)
    y <- runif(1, min=0, max=M)
    if (y < f(x, ...)) found=TRUE
  }
  return(x)
  #  arsample.unif(f, M, ...)
}

#' Symmetric Cayley distribution for angular data
#' 
#' The symmetric Cayley distribution has a density of the form \deqn{C_\mathrm{C}(r |\kappa)=\frac{1}{\sqrt{\pi}} \frac{\Gamma(\kappa+2)}{\Gamma(\kappa+1/2)}2^{-(\kappa+1)}(1+\cos r)^\kappa(1-\cos r)}.
#' It was orignally given in the material sciences literature by Schaben 1997 and called the de la Vall\'{e}e Poussin distribution but was more recently discussed and 
#' introduced in a more general manner by Leon 06.
#'
#' @param r Where the density is being evaluated
#' @param kappa The concentration paramter, taken to be zero
#' @param Haar logical, if density is evaluated with respect to Haar measure or Lebesgue
#' @return value of Cayley distribution with concentration \eqn{\kappa} evaluated at r
#' @seealso \code{\link{rcayley}},\code{\link{dfisher}},\code{\link{dhaar}}
#' @cite Schaeben97 leon06

dcayley <-function(r,kappa=1,Haar=F){
  den<-.5*gamma(kappa+2)/(sqrt(pi)*2^kappa*gamma(kappa+.5))*(1+cos(r))^kappa*(1-cos(r))
  
  if(Haar)
    return(den/(1-cos(r)))
  else
    return(den)
}

#' von Mises-Fisher distribution for angular data
#' 
#' The symmetric matrix fisher distribution has the density\deqn{C_\mathrm{{F}}(r|\kappa)=\frac{1}{2\pi[\mathrm{I_0}(2\kappa)-\mathrm{I_1}(2\kappa)]}e^{2\kappa\cos(r)}[1-\cos(r)]} 
#' where \eqn{\mathrm{I_p}(\cdot)} denotes the Bessel function of order \eqn{p} defined as  \eqn{\mathrm{I_p}(\kappa)=\frac{1}{2\pi}\int_{-\pi}^{\pi}\cos(pr)e^{\kappa\cos r}dr}.
#' This function allows the user to evaluate the function \eqn{C_\mathrm{{F}}(r|\kappa)} at \eqn{r} with \eqn{\kappa} provided by the user.
#'
#' @param r Where the density is being evaluated
#' @param kappa The concentration paramter, taken to be zero
#' @param Haar logical, if density is evaluated with respect to Haar measure or Lebesgue
#' @return value of Fisher matrix distribution with concentration \eqn{\kappa} evaluated at r

dfisher <-function(r,kappa=1,Haar=F){
  den<-exp(2*kappa*cos(r))*(1-cos(r))/(2*pi*(besselI(2*kappa,0)-besselI(2*kappa,1)))
  
  if(Haar)
    return(den/(1-cos(r)))
  
  else
    return(den)
}


#' Evaluate the uniform distribution on the circle at \eqn{r}
#' 
#' The uniform distribution on the sphere is also know as the Haar measure and has the density function \deqn{C_U(r)=\frac{1-cos(r)}{2\pi}}
#'
#' @param r Where the density is being evaluated
#' @return the probability density evaluated at r

dhaar <- function(r) return((1-cos(r))/(2*pi))

#'Density function for circular von Mises distribution
#'
#' The circular von Mises-based distribution has the density \deqn{C_\mathrm{M}(r|\kappa)=\frac{1}{2\pi \mathrm{I_0}(\kappa)}e^{\kappa\cos(r)}}.  This function allows the use to 
#' evaluate \eqn{C_\mathrm{M}(r|\kappa)} at angle \eqn{r} given a concentration parameter \eqn{\kappa}.
#'
#' @param r value at which to evaluate the distribution function
#' @param kappa concentration paramter
#' @param Haar logical, if density is evaluated with respect to Haar measure or Lebesgue
#' @return value of circular-von Mises distribution with concentration \eqn{\kappa} evaluated at r
#' @seealso \code{\link{rvmises}}, \code{\link{dfisher}},\code{\link{dhaar}},\code{\link{dcayley}}

dvmises<-function(r, kappa=1, Haar=F){
  den<-1/(2*pi*besselI(kappa,0))*exp(kappa*cos(r))
  
  if(Haar)
    return(den/(1-cos(r)))
  
  else
    return(den)
}

#' A function that will take in a Euler angle and return a rotation matrix in vector format
#' 
#' @param eur numeric Euler angle representation of an element in SO(3)
#' @return numeric \eqn{9\times 1} vector of a matrix in SO(3)
#' @seealso \code{\link{is.SO3}} can be used to check the output of this function
#' @export
#' @examples
#' eaExample<-c(pi/2,3*pi/4,0)
#' SO3Dat<-EAtoSO3(eaExample)
#' is.SO3(SO3Dat)

EAtoSO3<-function(eur){
  
  S<-matrix(NA,3,3)
  S[1,1]<-cos(eur[1])*cos(eur[3])-sin(eur[1])*sin(eur[3])*cos(eur[2])
  S[1,2]<-sin(eur[1])*cos(eur[3])+cos(eur[1])*sin(eur[3])*cos(eur[2])
  S[1,3]<-sin(eur[3])*sin(eur[2])
  
  S[2,1]<- -cos(eur[1])*sin(eur[3])-sin(eur[1])*cos(eur[3])*cos(eur[2])
  S[2,2]<- -sin(eur[1])*sin(eur[3])+cos(eur[1])*cos(eur[3])*cos(eur[2])
  S[2,3]<-cos(eur[3])*sin(eur[2])
  
  S[3,1]<-sin(eur[1])*sin(eur[2])
  S[3,2]<--cos(eur[1])*sin(eur[2])
  S[3,3]<-cos(eur[2])
  return(as.vector(S))
}

#' Directional vector to skew-symmetric Matrix
#' 
#' @author Heike Hofmann
#' @param U three dimensional vector indicating rotational fix-axis
#' @return skew-symmetric matrix 

eskew <- function(U) {
  U <- U/sqrt(sum(U^2))
  u <- U[1]
  v <- U[2]
  w <- U[3]
  
  res <- matrix((-1) *c(0,-w,v,w,0,-u,-v,u,0), ncol=3)
  return(res)
}

#' A novel approach to visualizing random rotations.
#' 
#' This function produces a three-dimensional globe onto which the on column of the provided sample is drawn.  The data are centered around a provided
#' matrix and the user can choose to display this center or not.  Based on \code{ggplot2} package by \cite{wickham09}.
#' 
#' @param Rs the sample of n random rotations
#' @param center point about which to center the observations
#' @param column integer 1 to 3 indicating which column to display
#' @param show.estimates rather to display the four estimates of the principal direction or not
#' @param ... Additional arguments passed to ggplot2
#' @return  a ggplot2 object with the data dispalyed on a blank sphere
#' @cite wickham09
#' @export
#' @examples
#' r<-rvmises(20,1.0)
#' Rs<-genR(r)
#' eyeBall(Rs,center=mean(Rs),show.estimates=TRUE,shape=4)

eyeBall<-function(Rs,center=diag(1,3,3),column=1,show.estimates=FALSE,...){
  
  # construct helper grid lines for sphere
  
  theta <- seq(0,pi, by=pi/8)
  phi <- seq(0,2*pi, by=0.005)
  df <- data.frame(expand.grid(theta=theta, phi=phi))
  
  #qplot(theta,phi, geom="point", data=df) + coord_polar()
  
  x <- with(df, sin(theta)*cos(phi))
  y <- with(df, sin(theta)*sin(phi))
  z <- with(df, cos(theta))
  circles <- data.frame(cbind(x,y,z))
  circles$ID <- as.numeric(factor(df$theta))
  
  theta <- seq(0,pi, by=0.005)
  phi <- seq(0,2*pi, by=pi/8)
  df <- data.frame(expand.grid(theta=theta, phi=phi))
  
  x <- with(df, sin(theta)*cos(phi))
  y <- with(df, sin(theta)*sin(phi))
  z <- with(df, cos(theta))
  circles.2 <- data.frame(cbind(x,y,z))
  circles.2$ID <- as.numeric(factor(df$phi))+9
  
  circles <- rbind(circles, circles.2)
  
  rot <- angle_axis(c(1,-1,0), pi/8)
  
  pcircles <- data.frame(as.matrix(circles[,1:3]) %*% rot)
  
  # this is the coordinate system and should be fixed, no matter what column of the rotation matrices is shown
  
  base <- ggplot(aes(x=X1, y=X2), data=pcircles[order(pcircles$X3), ]) + coord_equal() +  opts(legend.position="none") +
    geom_point(aes(colour=X3), size=0.6) + scale_colour_continuous(low=I("white"), high=I("grey50")) + opts(panel.background=theme_blank(),
    panel.grid.minor=theme_blank(),
    panel.grid.major=theme_blank(),
    axis.title.x=theme_blank(),
    axis.title.y=theme_blank(),
    axis.text.x=theme_blank(),
    axis.text.y=theme_blank(),
    axis.ticks=theme_blank())
  
  if(column==1){
    cols<-1:3
    rot<-angle_axis(c(0,1,0), pi/2) %*% rot
    
  }else if (column==2){
    cols<-4:6
    rot<-angle_axis(c(1,0,0), -pi/2) %*% rot
    
  }else{
    cols<-7:9
  }
  
  obs<-data.frame(as.matrix(Rs[,cols]) %*% center %*% rot)
  
  if(show.estimates){
    
    GMean<-as.vector(mean(Rs,type='intrinsic'))
    GMed<-as.vector(median.SO3(Rs,type='intrinsic'))
    PMed<-as.vector(median(Rs))
    PMean<-as.vector(mean(Rs))
    ests<-rbind(PMean,GMean,GMed,PMed)
    
    EstsDot<-data.frame(as.matrix(ests[,cols]) %*% center %*% rot)
    EstsDot$Shape<-as.factor(2:5)
    EstsDot$names<-labels(EstsDot)[[1]]
    
    out<-base+geom_point(aes(x=X1, y=X2,colour=X3),data=obs,...)+geom_point(aes(x=X1, y=X2,shape=Shape),data=EstsDot,size=2,...)
    
    out<-out+geom_text(aes(x=X1,y=X2,label=names),data=EstsDot,hjust=0,vjust=0)  }
  else{
    out<-base+geom_point(aes(x=X1, y=X2,colour=X3),data=obs,...)
  }
  return(out)
}



#' Generate rotation matrix given misorientation angle, r
#'
#' A function that generates a random rotation in \eqn{SO(3)} following a Uniform-Axis random roation distribution with central direction S
#' The exact form of the UARS distribution depends upon the distribution of the roation r
#'
#' @param r The angle through which all three dimensions are rotated after the axis was picked uniformly on the unit sphere
#' @param S the principle direction
#' @return a \eqn{n\times 9} matrix in SO(3) with misorientation angle r and principal direction S
#' @export
#' @examples
#' r<-rvmises(20,0.01)
#' genR(r)

genR<-function(r, S=diag(1,3,3)){
  
  if(!is.SO3(S)){
    stop("The principal direction must be in SO(3).")
  }
  
  o<-matrix(NA,length(r),9)
  I<-diag(1,3,3)
  
  #Generate angles theta from a uniform distribution from 0 to pi
  theta<-acos(runif(length(r),-1,1))
  
  #Generate angles phi from a uniform distribution from 0 to 2*pi
  phi<-runif(length(r),-pi,pi)
  
  for(i in 1:length(r)){
    
    #Using theta and phi generate a point uniformly on the unit sphere
    u<-matrix(c(sin(theta[i])*cos(phi[i]),sin(theta[i])*sin(phi[i]),cos(theta[i])),nrow=3,ncol=1,byrow=T)
    
    
    #Put it all together to make a rotation matrix O		
    o[i,]<-as.vector(S%*%angle_axis(u,r[i]))
    
  }
  class(o)<-"SO3"
  return(o)
}


#' A function to determine if a given matrix is in \eqn{SO(3)} or not.
#' 
#' @param x numeric \eqn{n \times n} matrix or vector of length \eqn{n^2}
#' @return logical T if the matrix is in SO(3) and false otherwise
#' @export
#' @examples
#' is.SO3(diag(1,3,3))
#' is.SO3(1:9)
is.SO3<-function(x){
  
  x<-matrix(x,3,3)
  
  #Does it have determinant 1?
  if(round(det(x),digits=7)!=1){
    return(FALSE)
  }
  
  #Is its transpose (approximately) its inverse?
  if(round(sum(t(x)%*%x-diag(1,3)),digits=10)!=0){
    return(FALSE)
  }
  
  return(TRUE)
  
}


#' This fuction will compute the natural exponential of skew-symmetric matrix.
#' 
#' See \cite{moakher02}
#' 
#' @param A 3-dimensional skew-symmetric matrix, i.e., \eqn{\bm A=-\bm A^\top}
#' @return numeric matrix \eqn{e^{\bm A}}
#' @cite moakher02

exp.skew<-function(A){
  
  if(round(sum(A-t(A)),digits=7)!=0){
    stop("The input matrix must be skew symmetric.")
  }
  
  I<-diag(1,3,3)
  AtA<-t(A)%*%A
  a2<-.5*sum(diag(AtA))
  a<-sqrt(a2)
  
  if(a==0){
    return(I)
  }else{
    p1<-(sin(a)/a)*A
    p2<-((1-cos(a))/a^2)*A%*%A
    return(I+p1+p2)
  }
  
}


#' This fuction will compute the natural logarithm of a matrix in SO(n).  It uses the special case of the Taylor expansion for SO(n) matrices.
#' 
#' For details see \cite{moakher02}
#' 
#' @param R numeric matrix in \eqn{SO(n)}
#' @return mlog numeric matrix \eqn{\log(R)}
#' @cite moakher02

log.SO3<-function(R){
  
  if(!is.SO3(R)){
    stop("This the input matrix must be in SO(n).")
  }
  
  tR<-sum(diag(R))
  cost<-.5*(tR-1)
  if(abs(cost)>=1){
    return(diag(0,3,3))
  }else{
    theta<-acos(cost)
    c<-theta/(2*sin(theta))
    mlog<-c*(R-t(R))
    return(mlog)
  }
}

#' Compute the projected or intrinsic mean estimate of the central direction
#'
#' This function takes a sample of \eqn{3\times 3} rotations (in the form of a \eqn{n\times 9} matrix where n is the sample size) and returns the projected arithmetic mean denoted \eqn{\widehat{\bm S}_P} or
#' intrinsic mean \eqn{\widehat{\bm S}_G} according to the \code{type} option.
#' For a sample of \eqn{n} random rotations \eqn{\bm{R}_i\in SO(3)$, $i=1,2,\dots,n}, the mean-type estimator is defined as \deqn{\widehat{\bm{S}}=\argmin_{\bm{S}\in SO(3)}\sum_{i=1}^nd_D^2(\bm{R}_i,\bm{S})} where \eqn{\bar{\bm{R}}=\frac{1}{n}\sum_{i=1}^n\bm{R}_i} and the distance metric \eqn{d_D}
#' is the Riemannian or Euclidean.  For more on the projected mean see \cite{moakher02} and for the intrinsic mean see \cite{manton04}.
#'
#' @param Rs A sample of n \eqn{3\times 3} random rotations
#' @param type String indicating 'projeted' or 'intrinsic' type mean estimator
#' @param epsilon Stopping rule for the intrinsic method
#' @param maxIter The maximum number of iterations allowed before returning most recent estimate
#' @return projected or intrinsic mean of the sample
#' @seealso \code{\link{median.SO3}}
#' @cite moakher02, manton04
#' @export
#' @examples
#' r<-rvmises(20,0.01)
#' Rs<-genR(r)
#' mean(Rs)

mean.SO3<-function(Rs, type='projected',epsilon=1e-5,maxIter=2000){
  
  if(!all(apply(Rs,1,is.SO3)))
    warning("Atleast one of the given observations is not in SO(3).  Use result with caution.")
  
  if(type != 'projected' & type!='intrinsic')
    stop("Incorrect usage of type option.  Select from 'projected' or 'intrinsic'.")
  
  R<-project.SO3(matrix(colMeans(Rs),3,3))
  
  if(type=='intrinsic'){
    
    n<-nrow(Rs)
    d<-1
    iter<-0
    s<-matrix(0,3,3)
    
    while(d>=epsilon){
      
      R<-R%*%exp.skew(s)
      
      s<-matrix(colMeans(t(apply(Rs,1,tLogMat,S=R))),3,3)
      
      d<-norm(s,type="F")
      
      iter<-iter+1
      
      if(iter>=maxIter){
        warning(paste("A unique solution wasn't found after",iter,"iterations."))
        return(R)
      }
    }
    
  }
  
  return(R)
}


#' Compute the projected or intrinsic median estimate of the central direction
#'
#' The median-type estimators are defined as \deqn{\widetilde{\bm{S}}=\argmin_{\bm{S}\in SO(3)}\sum_{i=1}^nd_D(\bm{R}_i,\bm{S})}.  If the choice of distance metrid, \eqn{d_D}, is Riemannian then the estimator is called the intrinsic, and if the distance metric in Euclidean then it projected.
#' The algorithm used in the intrinsic case is discussed in \cite{hartley11} and the projected case was written by the authors.
#' 
#' @param Rs the sample \eqn{n \times 9} matrix with rows corresponding to observations
#' @param type String indicating 'projeted' or 'intrinsic' type mean estimator
#' @param epsilon the stopping rule for the iterative algorithm 
#' @param maxIter integer, the maximum number of iterations allowed
#' @return S the element in SO(3) minimizing  the sum of first order Euclidean or Riemannian distances for sample Rs
#' @seealso \code{\link{mean.SO3}}
#' @cite hartley11
#' @export
#' @examples
#' r<-rcayley(50,1)
#' Rs<-genR(r)
#' median(Rs)

median.SO3<-function(Rs, type='projected',epsilon=1e-5,maxIter=2000){
  
  if(!all(apply(Rs,1,is.SO3)))
    warning("Atleast one of the given observations is not in SO(3).  Use result with caution.")
  
  if(type != 'projected' & type!='intrinsic')
    stop("Incorrect usage of type option.  Select from 'projected' or 'intrinsic'.")
  
  S<-mean(Rs)
  d<-1
  iter<-1
  delta<-matrix(0,3,3)
  
  while(d>=epsilon){
    
    if(type=='projected'){
      vn<-apply(Rs,1,vecNorm,type="F",S=S)
      
      delta<-matrix(colSums(Rs/vn)/sum(1/vn),3,3)
      
      Snew<-project.SO3(delta)
      
      d<-norm(Snew-S,type="F")
      S<-Snew
      
    }else if(type=='intrinsic'){
      
      S<-exp.skew(delta)%*%S
      
      v<-t(apply(Rs,1,tLogMat,S=S))
      vn<-apply(v,1,vecNorm,S=diag(0,3,3),type='F')
      
      delta<-matrix(colSums(v/vn)/sum(1/vn),3,3)
      d<-norm(delta,type="F")
      
    }
    
    iter<-iter+1
    
    if(iter>=maxIter){
      warning(paste("Unique solution wasn't found after ", iter, " iterations."))  
      return(S)
    }
  }
  return(S)
}

#' The projection of an arbitrary \eqn{3\times 3} matrix into \eqn{SO(3)}
#' 
#' This function uses the process given in Moakher 2002  to project an arbitrary \eqn{3\times 3} matrix into \eqn{SO(3)}.
#' @param M \eqn{3\times 3} matrix to project
#' @return projection of \eqn{\bm M} into \eqn{SO(3)}
#' @seealso \code{\link{mean.SO3}}, \code{\link{median.SO3}}
#' @export
#' @examples
#' M<-matrix(rnorm(9),3,3)
#' project.SO3(M)

project.SO3<-function(M){
  
  d<-svd(t(M)%*%M)
  u<-d$u
  
  d1<-1/sqrt(d$d[1])
  d2<-1/sqrt(d$d[2])
  d3<-sign(det(M))/sqrt(d$d[3])
  
  R<-M%*%u%*%diag(x=c(d1,d2,d3),3,3)%*%t(u)
  return(R)
}

#' A function to translate from unit quaternion representation to \eqn{SO(3)} representation
#' of a rotation matrix
#' 
#' @param q numeric unit vector, i.e. \eqn{q^\top q=1}, representing an element in SO(3)
#' @return vector representation of a rotation matrix in SO(3)
#' @seealso \code{\link{is.SO3}} can be used to check the return vector
#' @export
#' @examples
#' is.SO3(QtoSO3(c(1/sqrt(2),0,0,1/sqrt(2))))

QtoSO3<-function(q){
  
  if(round(t(q)%*%q,digits=10)!=1){
    stop("Input must have unit length.")
  }
  
  a<-q[1]
  b<-q[2]
  c<-q[3]
  d<-q[4]
  S<-matrix(NA,3,3)
  S[1,1]<-a^2+b^2-c^2-d^2
  S[1,2]<-2*b*c-2*a*d
  S[1,3]<-2*b*d+2*a*c
  S[2,1]<-2*b*c+2*a*d
  S[2,2]<-a^2-b^2+c^2-d^2
  S[2,3]<-2*c*d-2*a*b
  S[3,1]<-2*b*d-2*a*c
  S[3,2]<-2*c*d+2*a*b
  S[3,3]<-a^2-b^2-c^2+d^2
  
  return(as.vector(S))
}


#' Sample of size n from target density f
#' 
#' @author Heike Hofmann
#' @param n number of sample wanted
#' @param f target density
#' @param g sampling distribution
#' @param M maximum number in uniform proposal density
#' @param ... additional arguments sent to arsample
#' @return a vector of size n of observations from target density
#' @examples
#' # sample from haar distribution
#' x <- rar(10000, haar, runif, 1/pi, min=-pi, max=pi)
#' 
#' kappa=0.5
#' M <- max(fisher(seq(-pi, pi, length=1000), kappa))
#' x.fisher <- rar(10000, fisher, runif, M, min=-pi, max=pi, kappa=kappa)

rar <- function(n, f,g, M, ...) {
  res <- vector("numeric", length=n)
  for (i in 1:n) res[i] <- arsample(f,g,M, ...)
  return(res)
}

#' Simulate misorientation angles from Cayley distribtuion
#' 
#' This function allows the user to simulate \eqn{n} misorientation angles from the Cayley distribution symmetric about 0 on interval \eqn{(-\pi,\pi]}.  The relationship between Cayley and Beta distribution is used.
#' The symmetric Cayley distribution has a density of the form \deqn{C_\mathrm{C}(r |\kappa)=\frac{1}{\sqrt{\pi}} \frac{\Gamma(\kappa+2)}{\Gamma(\kappa+1/2)}2^{-(\kappa+1)}(1+\cos r)^\kappa(1-\cos r)}.
#' It was orignally given in the material sciences literature by Schaben 1997 and called the de la Vall\'{e}e Poussin distribution but was more recently discussed and 
#' introduced in a more general manner by Leon 06.
#' 
#' @param n sample size
#' @param kappa The concentration paramter
#' @return vector of n observations from Cayley(kappa) distribution
#' @cite Schaeben97 leon06
#' @export
#' @examples
#' r<-rcayley(20,0.01)

rcayley<-function(n,kappa=1){
  bet<-rbeta(n,kappa+0.5,3/2)
  theta<-acos(2*bet-1)*(1-2*rbinom(n,1,.5))
  return(theta)
}

#' Simulate a data set of size \eqn{n} from the matrix Fisher angular distribution
#' 
#' The symmetric matrix fisher distribution has the density\deqn{C_\mathrm{{F}}(r|\kappa)=\frac{1}{2\pi[\mathrm{I_0}(2\kappa)-\mathrm{I_1}(2\kappa)]}e^{2\kappa\cos(r)}[1-\cos(r)]} 
#' where \eqn{\mathrm{I_p}(\cdot)} denotes the Bessel function of order \eqn{p} defined as  \eqn{\mathrm{I_p}(\kappa)=\frac{1}{2\pi}\int_{-\pi}^{\pi}\cos(pr)e^{\kappa\cos r}dr}.
#' This function allows for simulation of \eqn{n} random deviates with density \eqn{C_\mathrm{{F}}(r|\kappa)} and \eqn{\kappa} provided by the user.
#' 
#' @param n sample size
#' @param kappa the concentration parameter
#' @return a sample of size \eqn{n} from the matrix Fisher distribution with concentration \eqn{\kappa}
#' @seealso \code{\link{dfisher}},\code{\link{rvmises}},\code{\link{rcayley}}


rfisher<-function(n,kappa=1){
  M<-max(dfisher(seq(-pi, pi, length=1000), kappa))
  return(rar(n,dfisher, runif, M, min=-pi, max=pi, kappa=kappa))
}

#' Riemannian Distance Between Two Random Rotations
#' 
#' This function will calculate the riemannian distance between an estimate of the central direction (in matrix or vector form) and the central direction.  By default the central direction
#' is taken to be the identity matrix, but any matrix in SO(3) will work.  It calls the matrix log and matrix exponential functions also given here.
#'
#' @param R The estimate of the central direction
#' @param S The true central direction
#' @return S3 \code{riedist} object; a number between 0 and pi that is the shortest geodesic curve connecting two matrices, i.e., the Riemannian distance
#' @export
#' @examples
#' r<-rvmises(20,0.01)
#' Rs<-genR(r)
#' Sp<-mean(Rs)
#' riedist(Sp,diag(1,3,3))

riedist<-function(R,S=diag(1,3,3)){
  R<-matrix(R,3,3)
  lRtS<-log.SO3(R%*%t(S))
  no<-norm(lRtS,type='F')
  return(no/sqrt(2))
}

#' Generate a vector of angles(r) from the von Mises Circular distribution
#'
#' The circular von Mises-based distribution has the density \deqn{C_\mathrm{M}(r|\kappa)=\frac{1}{2\pi \mathrm{I_0}(\kappa)}e^{\kappa\cos(r)}}.  This function allows the use to 
#' simulate \eqn{n} random deviates from \eqn{C_\mathrm{M}(r|\kappa)} given a concentration parameter \eqn{\kappa}.
#'
#' @param kappa The concentration parameter of the distribution
#' @param n The number of angles desired
#' @return S3 \code{rvmises} object; a vector of n angles following the von Mises Circular distribution with concentration kappa and mean/mode 0
#' @export
#' @examples
#' r<-rvmises(20,0.01)

rvmises<-function(n,kappa=1){
  u<-runif(3,0,1)
  a<-1+sqrt(1+4*kappa^2)
  b<-(a-sqrt(2*a))/(2*kappa)
  r<-(1+b^2)/(2*b)
  theta<-rep(10,n)
  
  for(i in 1:n){
    
    while(theta[i]==10){
      #Step  1
      u<-runif(3,0,1)
      z<-cos(pi*u[1])
      f<-(1+r*z)/(r+z)
      c<-kappa*(r-f)
      
      #Step 2
      u<-runif(3,0,1)
      if((c*(2-c)-u[2])>0){
        theta[i]=sign(u[3]-.5)*acos(f)
      }
      
      #Step 3
      else{
        if(log(c/u[2])+1-c<0){
          u<-runif(3,0,1)
        }
        
        else{
          u<-runif(3,0,1)
          theta[i]=sign(u[3]-.5)*acos(f)
        }
      }
    }
  }
  return(theta)
}


#' Compute the sum of the \eqn{p^{\text{th}}} order distances between Rs and S
#' 
#' @param Rs numeric matrix with sample size n rows and m columns
#' @param S the matrix to compute the sum of distances between each row of Rs with
#' @param p the order of the distances to compute
#' @return list of size two
#'  \item{Rieman}{the sum of \eqn{p^{\text{th}}} order Riemannian distances}
#'  \item{Euclid}{the sum of \eqn{p^{\text{th}}} order Euclidean distances}
#' @export
#' @examples
#' r<-rvmises(20,0.01)
#' Rs<-genR(r)
#' Sp<-mean(Rs)
#' SumDist(Rs,S=Sp,2)

SumDist<-function(Rs,S=diag(1,3,3),p){
  
  dist<-0
  dist2<-0
  n<-nrow(Rs)
  
  dR<-sum(apply(Rs,1,riedist,S=S)^p)
  
  dE<-sum(apply(Rs,1,vecNorm,type="F",S=S)^p)
  
  return(list(Rieman=dR,Euclid=dE))
}


tLogMat<-function(x,S){
  tra<-log.SO3(t(S)%*%matrix(x,3,3))
  return(as.vector(tra))
}


vecNorm<-function(x,S,...){
  n<-sqrt(length(x))
  cenX<-x-as.vector(S)
  return(norm(matrix(cenX,n,n),...))
}



qu <- function(Rs) {
  # represent rotation as quaternion
  theta <- eangle(Rs)
  u <- eaxis(Rs)
  x <- c(cos(theta/2), sin(theta/2)*u)
  names(x) <- c("s","i","j","k")
  x
}

euler <- function(rot) {
  # rotations to Euler angles
  if (is.matrix(rot)) rot <- as.vector(rot)
  alpha <- acos(-rot[8]/sqrt(1-rot[9]^2))
  beta <- acos(rot[9])
  gamma <- acos(rot[6]/sqrt(1-rot[9]^2))
  return(cbind(alpha, beta, gamma))
}
