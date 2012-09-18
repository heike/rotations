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

#' Method for creating a rotation using the angle axis representation
#'
#' Angle-axis representation based on the Rodrigues formula.
#'
#' @export
#' @param U three-dimensional vector describing the fix axis of the rotation
#' @param theta angle between -pi and pi
#' @return Used in \code{\link{eyeBall}} to orient the data properly

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


#' Accept/reject algorithm random sampling from angle distributions
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

#' Accept/reject algorithm random sampling from angle distributions using uniform envelop
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

#' Convert anything into EA class
#' 
#' @param x can be anything
#' @return x with class "EA"
#' @export

as.EA<-function(x){
  class(x)<-"EA"
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

#' Convert anything into SO3 class
#' 
#' @param x can be anything
#' @return x with class "SO3"
#' @export

as.SO3<-function(x){
  class(x)<-"SO3"
  return(x)
}


#' Find confidence interval radius for central orientation estimate
#' 
#' Given a sample or rotations and estimator for the central orientation, a confidence interval radius will
#' be estimated based on a bootstrap procedure.  The level of confidence, bootstrap iterations and bootstrap
#' sample sizes are all options for the user to set.  The procedure is as follows:
#'           0) Estimate Shat from (R_1,...,R_n) in Rs
#'           1) Sample R_1,...,R_m from Rs with replacement
#'           2) Estimate ShatStar from bootstrap sample
#'           3) Compute That=d_r(Shat,ShatStar)
#'           4) Repeat steps 1-3 B times
#'           5) Report q% percentile of That to be the radius of the confidence 'cone'
#'  @param Rs sample of size n rotations in matrix format
#'  @param main The method of central orientation estimation; e.g. mean or median, other methods can be used, but will result in a warning
#'  @param B bootstrap iterations
#'  @param m bootstrap sample size
#'  @param alpha level of confidence
#'  @param ... additional arguments passed to 'fun' such as 'type' and 'epsilon'
#'  @return the radius of the confidence cone
#'  @references bingham10
#'  @export
#'  @examples
#'  Rs<-genR(rvmises(20))
#'  CIradius(Rs,mean,type='projected')  ## calls CIradius.SO3
#'  Qs<-genR(rvmises(20),space='Q4')
#'  CIradius(Qs,mean,type='projected')  ## calls CIradius.Q4


CIradius <- function(Rs,main=mean,B=1000,m=nrow(Rs),alpha=0.95,...)
{
  UseMethod( "CIradius" )
}


#' @return \code{NULL}
#'
#' @rdname CIradius
#' @method CIradius SO3
#' @S3method CIradius SO3

CIradius.SO3 <- function(Rs,main=mean,B=1000,m=nrow(Rs),alpha=0.95,...){
  
  # This is more conservative than Bingham's method since d_r > max abs angle always
  args <- as.list(match.call())[-1]
  fname <- as.character(args$main)
  
#  if (!(fname %in% c("mean", "median")))
#    stop(sprintf("'%s' is not valid function for estimating main direction. Use 'mean' or 'median'.", fname))
  n<-nrow(Rs)
  Shat<-main(Rs,...)
  
  That <- rep(NA, B)
  for (i in 1:B) {  
    samp<-sample(1:n,m,replace=T)
    ShatStar<-main(as.SO3(Rs[samp,]),...)
    That[i] <- dist.SO3(Shat,ShatStar,method='intrinsic',p=1)
  }
  
  return(quantile(That,alpha))  
}

#' @return \code{NULL}
#'
#' @rdname CIradius
#' @method CIradius Q4
#' @S3method CIradius Q4

CIradius.Q4<-function(Qs,main=mean,B=1000,m=nrow(Qs),alpha=0.95,...){  
  Rs<-as.SO3(t(apply(Qs,1,SO3.Q4)))
#  args <- as.list(match.call())[-1]
#  fname <- as.character(args$main)
  
  return(CIradius(Rs,main,B,m,alpha,...))
}

#' @return \code{NULL}
#'
#' @rdname CIradius
#' @method CIradius EA
#' @S3method CIradius EA

CIradius.EA<-function(EAs,main=mean,B=1000,m=nrow(EAs),alpha=0.95,...){  
  Rs<-as.SO3(t(apply(EAs,1,SO3.EA)))
  #  args <- as.list(match.call())[-1]
  #  fname <- as.character(args$main)
  
  return(CIradius(Rs,main,B,m,alpha,...))
}

#' Symmetric Cayley distribution for angular data
#'
#' The symmetric Cayley distribution has a density of the form \deqn{C_\mathrm{C}(r |\kappa)=\frac{1}{\sqrt{\pi}} \frac{\Gamma(\kappa+2)}{\Gamma(\kappa+1/2)}2^{-(\kappa+1)}(1+\cos r)^\kappa(1-\cos r)}.
#' It was orignally given in the material sciences literature by Schaben 1997 and called the de la Vallee Poussin distribution but was more recently discussed and
#' introduced in a more general manner by Leon 06.
#'
#' @param r Where the density is being evaluated
#' @param kappa The concentration paramter, taken to be zero
#' @param Haar logical, if density is evaluated with respect to Haar measure or Lebesgue
#' @return value of Cayley distribution with concentration \eqn{\kappa} evaluated at r
#' @seealso \code{\link{rcayley}},\code{\link{dfisher}},\code{\link{dhaar}},\code{\link{dvmises}}
#' @export
#' @cite Schaeben97 leon06

dcayley <- function(r, kappa = 1, Haar = T) {
  den <- 0.5 * gamma(kappa + 2)/(sqrt(pi) * 2^kappa * gamma(kappa + 0.5)) * (1 + cos(r))^kappa * (1 - cos(r))
  
  if (Haar) 
    return(den/(1 - cos(r))) else return(den)
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
#' @seealso \code{\link{rfisher}}, \code{\link{dhaar}},\code{\link{dvmises}},\code{\link{dcayley}}
#' @export

dfisher <- function(r, kappa = 1, Haar = T) {
  den <- exp(2 * kappa * cos(r)) * (1 - cos(r))/(2 * pi * (besselI(2 * kappa, 0) - besselI(2 * kappa, 1)))
  
  if (Haar) {
    return(den/(1 - cos(r)))
  } else {
    return(den)
  }
}


#' Evaluate the uniform distribution on the circle at \eqn{r}
#'
#' The uniform distribution on the sphere is also know as the Haar measure and has the density function \deqn{C_U(r)=\frac{1-cos(r)}{2\pi}}
#'
#' @param r Where the density is being evaluated
#' @return the probability density evaluated at r
#' @seealso \code{\link{rhaar}}, \code{\link{dfisher}},\code{\link{dvmises}},\code{\link{dcayley}}
#' @export

dhaar <- function(r) return((1 - cos(r))/(2 * pi))

#' Distance Between Two Rotations
#'
#' This function will calculate the intrinsic (Riemannian) or projected (Euclidean) distance between two rotations.  If only one rotation is specified
#' the other will be set to the identity and the distance between the two is returned.
#'
#' @param R1 rotation in SO3 representation
#' @param R2 rotation in SO3 representation
#' @param method calculate intrinsic or projected distance
#' @param p the power of the respective distance
#' @return the pth power of the intrinsic or projected distance
#' @export
#' @examples
#' r<-rvmises(20,0.01)
#' Rs<-genR(r)
#' Sp<-mean(Rs)
#' dist.SO3(Sp)

dist.SO3 <- function(R1, R2=id.SO3, method='projected' , p=1) {
  
  R1<-matrix(R1,3,3)
  R2<-matrix(R2,3,3)
  
  if(method=='projected'){
    
    so3dist<-norm(R1-R2,type='F')^p
    
  }else if(method=='intrinsic'){
    
    so3dist<-eangle(t(R1)%*%R2)^p
    
  }else{
    stop("Incorrect usage of method argument.  Please choose intrinsic or projected.")
  }
  
  return(so3dist)
  
}

#' Distance Between Two Rotations
#'
#' This function will calculate the intrinsic (Riemannian) or projected (Euclidean) distance between two rotations.  If only one rotation is specified
#' the other will be set to the identity and the distance between the two is returned.
#'
#' @param Q1 First rotation in quaternion form
#' @param Q2 Second rotation in quaternion form, identity by default
#' @param method calculate intrinsic or projected distance
#' @param p the power of the respective distance
#' @return the pth power of the intrinsic or projected distance between Q1 and Q2
#' @export
#' @examples
#' r<-rvmises(20,0.01)
#' Qs<-genR(r,space="Q4")
#' Qp<-mean(Qs)
#' dist(Qp)

dist.Q4 <- function(Q1, Q2=id.Q4 ,method='projected', p=1) {
  
  
  if(method=='intrinsic'){
    
    cp <- sum(Q1*Q2)
    q4dist<-acos(2*cp*cp-1)^p
    
  }else if(method=='projected'){
    
    R1<-SO3.Q4(Q1)
    R2<-SO3.Q4(Q2)
    q4dist<-norm(R1-R2,type='F')^p
    
  }else{
    stop("Incorrect usage of method argument.  Please choose intrinsic or projected.")
  }
  
  return(q4dist)
}

#' Distance Between Two Rotations
#'
#' This function will calculate the intrinsic (Riemannian) or projected (Euclidean) distance between two rotations.  If only one rotation is specified
#' the other will be set to the identity and the distance between the two is returned.
#'
#' @param s The estimate of the central direction
#' @param method calculate intrinsic or projected distance
#' @param p the power of the respective distance
#' @return the pth power of the projected or intrinsic distance between Q1 and Q2
#' @export
#' @examples
#' r<-rvmises(20,0.01)
#' EAs<-genR(r,space="EA")
#' EAp<-mean(EAs)
#' dist(EAp)

dist.EA <- function(EA1, EA2=id.EA ,method='projected', p=1) {
  
  R1<-SO3.EA(EA1)
  
  R2<-as.SO3(matrix(SO3.EA(EA2),3,3))
  
  EAdist<-dist.SO3(R1,R2,method,p)
  
  return(EAdist)
}


#'Density function for circular von Mises distribution
#'
#' The circular von Mises-based distribution has the density \deqn{C_\mathrm{M}(r|\kappa)=\frac{1}{2\pi \mathrm{I_0}(\kappa)}e^{\kappa\cos(r)}}.  This function allows the use to
#' evaluate \eqn{C_\mathrm{M}(r|\kappa)} at angle \eqn{r} given a concentration parameter \eqn{\kappa}.
#'
#' @param r value at which to evaluate the distribution function
#' @param kappa concentration paramter
#' @param Haar logical, if density is evaluated with respect to Haar measure or Lebesgue
#' @export
#' @return value of circular-von Mises distribution with concentration \eqn{\kappa} evaluated at r
#' @seealso \code{\link{rvmises}}, \code{\link{dfisher}},\code{\link{dhaar}},\code{\link{dcayley}}

dvmises <- function(r, kappa = 1, Haar = T) {
  den <- 1/(2 * pi * besselI(kappa, 0)) * exp(kappa * cos(r))
  
  if (Haar) {
    return(den/(1 - cos(r)))
  } else {
    return(den)
  }
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
  
  #For now deal with singularity this was
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

#' Find the angle of rotation R
#' 
#' This is a test
#' This function will find the angle of rotation matrix R.  Used mostly to reparameterize between
#' quaternions and rotations.
#' @param R 3-by-3 matrix in SO3 whos angle of rotation is needed
#' @return angle of rotation
#' @seealso \code{\link{eaxis}}

eangle<-function(Rs){
  # for R in SO(3):
  #  1 + 2 cos(theta) = tr(R)
  
  tr<-Rs[1]+Rs[5]+Rs[9]
  return(acos((tr-1)/2))
}

#' Find the axis of rotation R
#' 
#' This function will find the axis of rotation matrix R.  The simple calculation is based on Rodrigues formula
#' and noticing that R - t(R) can be simplified greatly.  Used mostly to reparameterize between
#' quaternions and rotations.
#' @param R 3-by-3 matrix in SO3 whos axis is needed
#' @return three dimensional axis of length one
#' @seealso \code{\link{eangle}}

eaxis<-function(R){
  # based on Rodrigues formula: R - t(R)
  
  R<-matrix(R,3,3)
  
  X <- R - t(R)
  u <- rev(X[upper.tri(X)])*c(-1,1,-1)
  
  return(u/sqrt(sum(u^2))) # will be trouble, if R is symmetric, i.e. id,  .... 

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

#' A novel approach to visualizing random rotations.
#'
#' This function produces a three-dimensional globe onto which the on column of the provided sample is drawn.  The data are centered around a provided
#' matrix and the user can choose to display this center or not.  Based on \code{ggplot2} package by \cite{wickham09}.
#'
#' @param Rs the sample of n random rotations
#' @param center point about which to center the observations
#' @param column integer 1 to 3 indicating which column to display
#' @param show_estimates rather to display the four estimates of the principal direction or not
#' @param xlimits limits for the x-axis, appropriate range is [-1,1]
#' @param ylimits limits for the y-axis, appropriate range is [-1,1]
#' @param ... Additional arguments passed to ggplot2
#' @return  a ggplot2 object with the data dispalyed on a blank sphere
#' @cite wickham09
#' @export
#' @examples
#' r<-rvmises(20,1.0)
#' Rs<-genR(r)
#' eyeBall(Rs,center=mean(Rs),show_estimates=TRUE,shape=4)

eyeBall <- function(Rs, center = id.SO3, column = 1, show_estimates = FALSE, xlimits=c(-1,1),ylimits=c(-1,1), ...) {
  
  # construct helper grid lines for sphere
  
  theta <- seq(0, pi, by = pi/8)
  phi <- seq(0, 2 * pi, by = 0.005)
  df <- data.frame(expand.grid(theta = theta, phi = phi))
  
  # qplot(theta,phi, geom='point', data=df) + coord_polar()
  
  x <- with(df, sin(theta) * cos(phi))
  y <- with(df, sin(theta) * sin(phi))
  z <- with(df, cos(theta))
  circles <- data.frame(cbind(x, y, z))
  circles$ID <- as.numeric(factor(df$theta))
  
  theta <- seq(0, pi, by = 0.005)
  phi <- seq(0, 2 * pi, by = pi/8)
  df <- data.frame(expand.grid(theta = theta, phi = phi))
  
  x <- with(df, sin(theta) * cos(phi))
  y <- with(df, sin(theta) * sin(phi))
  z <- with(df, cos(theta))
  circles.2 <- data.frame(cbind(x, y, z))
  circles.2$ID <- as.numeric(factor(df$phi)) + 9
  
  circles <- rbind(circles, circles.2)
  
  rot <- SO3(c(1, -1, 0), pi/8)
  
  pcircles <- data.frame(as.matrix(circles[, 1:3]) %*% rot)
  
  # this is the coordinate system and should be fixed, no matter what column of the rotation matrices is
  # shown
  
  base <- ggplot(aes(x = X1, y = X2), data = pcircles[order(pcircles$X3), ]) + coord_equal() + opts(legend.position = "none") + 
    geom_point(aes(colour = X3), size = 0.6) + scale_colour_continuous(low = I("white"), high = I("grey50")) + 
    opts(panel.background = theme_blank(), panel.grid.minor = theme_blank(), panel.grid.major = theme_blank(), 
         axis.title.x = theme_blank(), axis.title.y = theme_blank(), axis.text.x = theme_blank(), axis.text.y = theme_blank(), 
         axis.ticks = theme_blank())+xlim(xlimits)+ylim(ylimits)
  
  if (column == 1) {
    cols <- 1:3
    rot <- SO3(c(0, 1, 0), pi/2) %*% rot
    
  } else if (column == 2) {
    cols <- 4:6
    rot <- SO3(c(1, 0, 0), -pi/2) %*% rot
    
  } else {
    cols <- 7:9
  }
  
  obs <- data.frame(as.matrix(Rs[, cols]) %*% center %*% rot)
  
  if (show_estimates) {
    
    GMean <- as.vector(mean(Rs, type = "intrinsic"))
    GMed <- as.vector(median.SO3(Rs, type = "intrinsic"))
    PMed <- as.vector(median(Rs))
    PMean <- as.vector(mean(Rs))
    ests <- rbind(PMean, GMean, GMed, PMed)
    
    EstsDot <- data.frame(as.matrix(ests[, cols]) %*% center %*% rot)
    EstsDot$Shape <- as.factor(2:5)
    EstsDot$names <- labels(EstsDot)[[1]]
    
    out <- base + geom_point(aes(x = X1, y = X2, colour = X3), data = obs, ...) + geom_point(aes(x = X1, 
                                                                                                 y = X2, shape = Shape), data = EstsDot, size = 2, ...)
    
    out <- out + geom_text(aes(x = X1, y = X2, label = names), data = EstsDot, hjust = 0, vjust = 0)
  } else {
    out <- base + geom_point(aes(x = X1, y = X2, colour = X3), data = obs, ...)
  }
  return(out)
}


eyeBallwCI <- function(Rs, center = id.SO3, column = 1, show_estimates = FALSE, xlimits=c(-1,1),ylimits=c(-1,1), ...) {
  
  # construct helper grid lines for sphere
  
  theta <- seq(0, pi, by = pi/8)
  phi <- seq(0, 2 * pi, by = 0.005)
  df <- data.frame(expand.grid(theta = theta, phi = phi))
  
  # qplot(theta,phi, geom='point', data=df) + coord_polar()
  
  x <- with(df, sin(theta) * cos(phi))
  y <- with(df, sin(theta) * sin(phi))
  z <- with(df, cos(theta))
  circles <- data.frame(cbind(x, y, z))
  circles$ID <- as.numeric(factor(df$theta))
  
  theta <- seq(0, pi, by = 0.005)
  phi <- seq(0, 2 * pi, by = pi/8)
  df <- data.frame(expand.grid(theta = theta, phi = phi))
  
  x <- with(df, sin(theta) * cos(phi))
  y <- with(df, sin(theta) * sin(phi))
  z <- with(df, cos(theta))
  circles.2 <- data.frame(cbind(x, y, z))
  circles.2$ID <- as.numeric(factor(df$phi)) + 9
  
  circles <- rbind(circles, circles.2)
  
  rot <- SO3(c(1, -1, 0), pi/8)
  
  pcircles <- data.frame(as.matrix(circles[, 1:3]) %*% rot)
  
  # this is the coordinate system and should be fixed, no matter what column of the rotation matrices is
  # shown
  
  base <- ggplot(aes(x = X1, y = X2), data = pcircles[order(pcircles$X3), ]) + coord_equal() + opts(legend.position = "none") + 
    geom_point(aes(colour = X3), size = 0.6) + scale_colour_continuous(low = I("white"), high = I("grey50")) + 
    opts(panel.background = theme_blank(), panel.grid.minor = theme_blank(), panel.grid.major = theme_blank(), 
         axis.title.x = theme_blank(), axis.title.y = theme_blank(), axis.text.x = theme_blank(), axis.text.y = theme_blank(), 
         axis.ticks = theme_blank())+xlim(xlimits)+ylim(ylimits)
  
  if (column == 1) {
    cols <- 1:3
    rot <- SO3(c(0, 1, 0), pi/2) %*% rot
    
  } else if (column == 2) {
    cols <- 4:6
    rot <- SO3(c(1, 0, 0), -pi/2) %*% rot
    
  } else {
    cols <- 7:9
  }
  
  obs <- data.frame(as.matrix(Rs[, cols]) %*% center %*% rot)
  
  if (show_estimates) {
    
    GMean <- as.vector(mean(Rs, type = "intrinsic"))
    GMeanRad<-CIradius(Rs)
    
    GMean.boot <- t(replicate(2000, SO3(c(runif(2,-1,1),0), GMeanRad),simplify="matrix"))
    
    GMean.sp<-data.frame(as.matrix(GMean.boot[,7:9]) %*% t(matrix(GMean,3,3)) %*% center %*% rot)
    
    GMed <- as.vector(median.SO3(Rs, type = "intrinsic"))
    PMed <- as.vector(median(Rs))
    PMean <- as.vector(mean(Rs))
    ests <- rbind(PMean, GMean, GMed, PMed)
    
    EstsDot <- data.frame(as.matrix(ests[, cols]) %*% center %*% rot)
    EstsDot$Shape <- as.factor(2:5)
    EstsDot$names <- labels(EstsDot)[[1]]
    
    out <- base + geom_point(aes(x = X1, y = X2, colour = X3), data = obs, ...) + geom_point(aes(x = X1, 
                                                                                                 y = X2, shape = Shape), data = EstsDot, size = 2, ...)
    
    out <- out + geom_text(aes(x = X1, y = X2, label = names), data = EstsDot, hjust = 0, vjust = 0)
    out <- out + geom_point(aes(x=X1,y=X2,colour= X3),data=GMean.sp)
    
  } else {
    out <- base + geom_point(aes(x = X1, y = X2, colour = X3), data = obs, ...)
  }
  return(out)
}



#' Generate rotation matrix given misorientation angle, r
#'
#' A function that generates a random rotation in \eqn{SO(3)} following a Uniform-Axis random roation distribution with central direction S
#' The exact form of the UARS distribution depends upon the distribution of the roation r
#'
#' @param r The angle through which all three dimensions are rotated after the axis was picked uniformly on the unit sphere
#' @param S The principle direction
#' @param space Indicates the desired representation: matrix in SO3, quaternion, or Euler angles 
#' @return a matrix where each row is a sample point in the desired space
#' @export
#' @examples
#' r<-rvmises(20,0.01)
#' genR(r)

genR <- function(r, S = diag(1, 3, 3), space='SO3') {
  
  if(!(space %in% c("SO3","Q4","EA")))
    stop("Incorrect space argument.  Options are: SO3, Q4 and EA. ")
  
  if (!is.SO3(S))
    stop("The principal direction must be in SO(3).")
  
  if(space=="SO3")
    o <- matrix(NA, length(r), 9)
  else if(space=="Q4")
    q <- matrix(NA, length(r), 4)
  else
    ea <- matrix(NA, length(r), 3)

  
  # Generate angles theta from a uniform distribution from 0 to pi
  theta <- acos(runif(length(r), -1, 1))
  
  # Generate angles phi from a uniform distribution from 0 to 2*pi
  phi <- runif(length(r), -pi, pi)
  
  for (i in 1:length(r)) {
    
    # Using theta and phi generate a point uniformly on the unit sphere
    u <- c(sin(theta[i]) * cos(phi[i]), sin(theta[i]) * sin(phi[i]), cos(theta[i]))
    
    if(space=="SO3"){
      
      o[i,] <- as.vector(S %*% SO3(u, r[i]))
      
    }else if(space=="Q4"){
      
      q[i,] <- c(cos(r[i]/2),sin(r[i]/2)*S%*%u)
      
    }else{
      
      ea[i,] <- EA.SO3(S %*% SO3(u, r[i]))
    }
  }
  if(space=="SO3"){
    class(o) <- "SO3"
    return(o)
  }else if (space=="Q4"){
    class(q)<-"Q4"
    return(q)
  }else{
    class(ea)<-"EA"
    return(ea)
  }
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
  
  # Does it have determinant 1?
  if (abs(det(x)- 1)>10e-10) {
    return(FALSE)
  }
  
  # Is its transpose (approximately) its inverse?
  if (sum(t(x) %*% x - diag(1, 3))>10e-10) {
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


#' This fuction will compute the natural logarithm of a matrix in SO(n).  It uses the special case of the Taylor expansion for SO(n) matrices.
#'
#' For details see \cite{moakher02}
#'
#' @param R numeric matrix in \eqn{SO(n)}
#' @return mlog numeric matrix \eqn{\log(R)}
#' @cite moakher02

log.SO3 <- function(R) {
  
  if (!is.SO3(R)) {
    stop("This the input matrix must be in SO(3).")
  }
  
  theta <- eangle(R)
  
  if (abs(cos(theta)) >= 1) {
    return(diag(0, 3, 3))
  } else {
    
    c <- theta/(2 * sin(theta))
    mlog <- c * (R - t(R))
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
#' @param Rs A \eqn{n\times 9} matrix where each row corresponds to a random rotation in matrix form
#' @param type String indicating 'projeted' or 'intrinsic' type mean estimator
#' @param epsilon Stopping rule for the intrinsic method
#' @param maxIter The maximum number of iterations allowed before returning most recent estimate
#' @return projected or intrinsic mean of the sample
#' @seealso \code{\link{median.SO3}}
#' @cite moakher02, manton04
#' @S3method mean SO3
#' @method mean SO3
#' @examples
#' r<-rvmises(20,0.01)
#' Rs<-genR(r)
#' mean(Rs)

mean.SO3 <- function(Rs, type = "projected", epsilon = 1e-05, maxIter = 2000) {
  if (!all(apply(Rs, 1, is.SO3))) 
    warning("At least one of the observations is not in SO(3).  Use result with caution.")
  
  if (!(type %in% c("projected", "intrinsic")))
    stop("type needs to be one of 'projected' or 'intrinsic'.")
  
  R <- project.SO3(matrix(colMeans(Rs), 3, 3))
  
  if (type == "intrinsic") {
    
    n <- nrow(Rs)
    d <- 1
    iter <- 0
    s <- matrix(0, 3, 3)
    
    while (d >= epsilon) {
      
      R <- R %*% exp.skew(s)
      
      s <- matrix(colMeans(t(apply(Rs, 1, tLogMat, S = R))), 3, 3)
      
      d <- norm(s, type = "F")
      
      iter <- iter + 1
      
      if (iter >= maxIter) {
        warning(paste("No convergence in ", iter, " iterations."))
        return(as.SO3(R))
      }
    }
    
  }
  class(R)<-"SO3"
  return(R)
}

#' Compute the projected or intrinsic mean estimate of the central direction
#'
#' This function takes a sample of n unit quaternions and approximates the mean rotation.  If the projected mean
#' is called for then the quaternions are turned reparameterized to matrices and mean.SO3 is called.  If the intrinsic
#' mean is called then according to \cite{gramkow01} a better approximation is achieved by taking average quaternion
#' and normalizing.  Our simulations don't match this claim.
#'
#' 
#' @param Qs A \eqn{n\times 4} matrix where each row corresponds to a random rotation in unit quaternion
#' @param type String indicating 'projeted' or 'intrinsic' type mean estimator
#' @param epsilon Stopping rule for the intrinsic method
#' @param maxIter The maximum number of iterations allowed before returning most recent estimate
#' @return projected or intrinsic mean of the sample
#' @seealso \code{\link{mean.SO3}}
#' @cite moakher02, manton04
#' @export
#' @examples
#' r<-rvmises(20,0.01)
#' Qs<-genR(r,space="Q4")
#' mean(Qs,type='intrinsic')

mean.Q4 <- function(Qs, type = "projected", epsilon = 1e-05, maxIter = 2000) {
  
  Rs<-t(apply(Qs,1,SO3.Q4))
  
  R<-mean.SO3(Rs,type,epsilon,maxIter)
    
  return(Q4.SO3(R))
    
}

#' Compute the projected or intrinsic mean estimate of the central direction
#'
#' This function takes a sample of \eqn{3\times 3} rotations (in the form of a \eqn{n\times 9} matrix where n is the sample size) and returns the projected arithmetic mean denoted \eqn{\widehat{\bm S}_P} or
#' intrinsic mean \eqn{\widehat{\bm S}_G} according to the \code{type} option.
#' For a sample of \eqn{n} random rotations \eqn{\bm{R}_i\in SO(3)$, $i=1,2,\dots,n}, the mean-type estimator is defined as \deqn{\widehat{\bm{S}}=\argmin_{\bm{S}\in SO(3)}\sum_{i=1}^nd_D^2(\bm{R}_i,\bm{S})} where \eqn{\bar{\bm{R}}=\frac{1}{n}\sum_{i=1}^n\bm{R}_i} and the distance metric \eqn{d_D}
#' is the Riemannian or Euclidean.  For more on the projected mean see \cite{moakher02} and for the intrinsic mean see \cite{manton04}.
#'
#' @param EAs A \eqn{n\times 3} matrix where each row corresponds to a random rotation in Euler angle form
#' @param type String indicating 'projeted' or 'intrinsic' type mean estimator
#' @param epsilon Stopping rule for the intrinsic method
#' @param maxIter The maximum number of iterations allowed before returning most recent estimate
#' @return projected or intrinsic mean of the sample
#' @seealso \code{\link{mean.SO3}}
#' @cite moakher02, manton04
#' @export
#' @examples
#' r<-rvmises(20,0.01)
#' EAs<-genR(r,space="R3")
#' mean.EA(EAs)

mean.EA <- function(EAs, type = "projected", epsilon = 1e-05, maxIter = 2000) {
  
  Rs<-as.SO3(t(apply(EAs,1,SO3.EA)))
  
  R<-mean(Rs,type,epsilon,maxIter)
  
  return(EA.SO3(R))
}


#' Compute the projected or intrinsic median estimate of the central direction
#'
#' The median-type estimators are defined as \deqn{\widetilde{\bm{S}}=\argmin_{\bm{S}\in SO(3)}\sum_{i=1}^nd_D(\bm{R}_i,\bm{S})}.  If the choice of distance metrid, \eqn{d_D}, is Riemannian then the estimator is called the intrinsic, and if the distance metric in Euclidean then it projected.
#' The algorithm used in the intrinsic case is discussed in \cite{hartley11} and the projected case was written by the authors.
#'
#' @param Rs A \eqn{n\times 9} matrix where each row corresponds to a random rotation in matrix form
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

median.SO3 <- function(Rs, type = "projected", epsilon = 1e-05, maxIter = 2000) {
  
  if (!all(apply(Rs, 1, is.SO3))) 
    warning("Atleast one of the given observations is not in SO(3).  Use result with caution.")
  
  if (type != "projected" & type != "intrinsic") 
    stop("Incorrect usage of type option.  Select from 'projected' or 'intrinsic'.")
  
  S <- mean(Rs)
  d <- 1
  iter <- 1
  delta <- matrix(0, 3, 3)
  
  while (d >= epsilon) {
    
    if (type == "projected") {
      vn <- apply(Rs, 1, vecNorm, type = "F", S = S)
      
      delta <- matrix(colSums(Rs/vn)/sum(1/vn), 3, 3)
      
      Snew <- project.SO3(delta)
      
      d <- norm(Snew - S, type = "F")
      S <- Snew
      
    } else if (type == "intrinsic") {
      
      S <- exp.skew(delta) %*% S
      
      v <- t(apply(Rs, 1, tLogMat, S = S))
      vn <- apply(v, 1, vecNorm, S = diag(0, 3, 3), type = "F")
      
      delta <- matrix(colSums(v/vn)/sum(1/vn), 3, 3)
      d <- norm(delta, type = "F")
      
    }
    
    iter <- iter + 1
    
    if (iter >= maxIter) {
      warning(paste("Unique solution wasn't found after ", iter, " iterations."))
  
      return(as.SO3(S))
    }
  }
  class(S)<-"SO3"
  return(S)
}

#' Compute the projected or intrinsic median estimate of the central direction
#'
#' The median-type estimators are defined as \deqn{\widetilde{\bm{S}}=\argmin_{\bm{S}\in SO(3)}\sum_{i=1}^nd_D(\bm{R}_i,\bm{S})}.  If the choice of distance metrid, \eqn{d_D}, is Riemannian then the estimator is called the intrinsic, and if the distance metric in Euclidean then it projected.
#' The algorithm used in the intrinsic case is discussed in \cite{hartley11} and the projected case was written by the authors.
#'
#' @param Qs A \eqn{n\times 4} matrix where each row corresponds to a random rotation in unit quaternion form
#' @param type String indicating 'projeted' or 'intrinsic' type mean estimator
#' @param epsilon the stopping rule for the iterative algorithm
#' @param maxIter integer, the maximum number of iterations allowed
#' @return S the unit quaternion minimizing  the sum of first order Euclidean or Riemannian distances for sample Rs
#' @seealso \code{\link{median.SO3}}
#' @cite hartley11
#' @export
#' @examples
#' r<-rcayley(50,1)
#' Qs<-genR(r,space="Q4")
#' median(Qs)

median.Q4 <- function(Qs, type = "projected", epsilon = 1e-05, maxIter = 2000) {

  Rs<-t(apply(Qs,1,SO3.Q4))
  
  R<-median.SO3(Rs,type,epsilon,maxIter)
  
  return(Q4.SO3(R))
}

#' Compute the projected or intrinsic median estimate of the central direction
#'
#' The median-type estimators are defined as \deqn{\widetilde{\bm{S}}=\argmin_{\bm{S}\in SO(3)}\sum_{i=1}^nd_D(\bm{R}_i,\bm{S})}.  If the choice of distance metrid, \eqn{d_D}, is Riemannian then the estimator is called the intrinsic, and if the distance metric in Euclidean then it projected.
#' The algorithm used in the intrinsic case is discussed in \cite{hartley11} and the projected case was written by the authors.
#'
#' @param EAs A \eqn{n\times 3} matrix where each row corresponds to a random rotation in Euler angle form
#' @param type String indicating 'projeted' or 'intrinsic' type mean estimator
#' @param epsilon the stopping rule for the iterative algorithm
#' @param maxIter integer, the maximum number of iterations allowed
#' @return S the vector of Euler angles minimizing  the sum of first order Euclidean or Riemannian distances for sample Rs
#' @seealso \code{\link{median.SO3}}
#' @cite hartley11
#' @export
#' @examples
#' r<-rcayley(50,1)
#' EA<-genR(r,space="R3")
#' median.R3(EA)

median.EA <- function(EAs, type = "projected", epsilon = 1e-05, maxIter = 2000) {
  
  Rs<-t(apply(EAs,1,SO3.EA))
  
  R<-median.SO3(Rs,type,epsilon,maxIter)
  
  return(EA.SO3(R))
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

project.SO3 <- function(M) {
  
  d <- svd(t(M) %*% M)
  u <- d$u
  
  d1 <- 1/sqrt(d$d[1])
  d2 <- 1/sqrt(d$d[2])
  d3 <- sign(det(M))/sqrt(d$d[3])
  
  R <- M %*% u %*% diag(x = c(d1, d2, d3), 3, 3) %*% t(u)
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
#' SO3Dat<-SO3(eaExample)
#' is.SO3(SO3Dat)

SO3.EA <- function(eur) {
  
  S <- matrix(NA, 3, 3)
  S[1, 1] <- cos(eur[1]) * cos(eur[3]) - sin(eur[1]) * sin(eur[3]) * cos(eur[2])
  S[1, 2] <- sin(eur[1]) * cos(eur[3]) + cos(eur[1]) * sin(eur[3]) * cos(eur[2])
  S[1, 3] <- sin(eur[3]) * sin(eur[2])
  
  S[2, 1] <- -cos(eur[1]) * sin(eur[3]) - sin(eur[1]) * cos(eur[3]) * cos(eur[2])
  S[2, 2] <- -sin(eur[1]) * sin(eur[3]) + cos(eur[1]) * cos(eur[3]) * cos(eur[2])
  S[2, 3] <- cos(eur[3]) * sin(eur[2])
  
  S[3, 1] <- sin(eur[1]) * sin(eur[2])
  S[3, 2] <- -cos(eur[1]) * sin(eur[2])
  S[3, 3] <- cos(eur[2])
  S <- as.SO3(as.vector(S))

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
  
  
  theta<-2*acos(q[1])
  
  if(theta==0){
    u <- rep(0,3)
  }else{
    u <- q[2:4]/sin(theta/2)
  }
  return(SO3(u, theta)) 
}

#' Reparameterize a rotation matrix as a unit quaternion
#' 
#' The rotation R is taken, the angle and axis of rotation are found and a unit quaternion is formed from the results.
#' @param R a rotation matrix in SO3
#' @return a unit quaternion of class Q4

Q4.SO3 <- function(R) {
  
  theta <- eangle(R)
  u <- eaxis(R)
  x <- as.Q4(c(cos(theta/2), sin(theta/2) * u))

  return(x)
}


#' Sample of size n from target density f
#'
#' @author Heike Hofmann
#' @param n number of sample wanted
#' @param f target density
#' @param M maximum number in uniform proposal density
#' @param ... additional arguments sent to arsample
#' @return a vector of size n of observations from target density

rar <- function(n, f, M, ...) {
  res <- vector("numeric", length = n)
  for (i in 1:n) res[i] <- arsample.unif(f, M, ...)
  return(res)
}


#' Simulate misorientation angles from Cayley distribtuion
#'
#' This function allows the user to simulate \eqn{n} misorientation angles from the Cayley distribution symmetric about 0 on interval \eqn{(-\pi,\pi]}.  The relationship between Cayley and Beta distribution is used.
#' The symmetric Cayley distribution has a density of the form \deqn{C_\mathrm{C}(r |\kappa)=\frac{1}{\sqrt{\pi}} \frac{\Gamma(\kappa+2)}{\Gamma(\kappa+1/2)}2^{-(\kappa+1)}(1+\cos r)^\kappa(1-\cos r)}.
#' It was orignally given in the material sciences literature by Schaben 1997 and called the de la Vallee Poussin distribution but was more recently discussed and
#' introduced in a more general manner by Leon 06.
#'
#' @param n sample size
#' @param kappa The concentration paramter
#' @param nu An alternative to kappa; circular variance
#' @return vector of n observations from Cayley(kappa) distribution
#' @cite Schaeben97 leon06
#' @seealso \code{\link{dcayley}},\code{\link{rvmises}},\code{\link{rcayley}},\code{\link{rhaar}}
#' @export
#' @examples
#' r<-rcayley(20,0.01)

rcayley <- function(n, kappa = 1, nu = NULL) {
  
  if(!is.null(nu))
    kappa <- cayley_kappa(nu)
      
  bet <- rbeta(n, kappa + 0.5, 3/2)
  theta <- acos(2 * bet - 1) * (1 - 2 * rbinom(n, 1, 0.5))
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
#' @param nu An alternative to kappa; circular variance
#' @return a sample of size \eqn{n} from the matrix Fisher distribution with concentration \eqn{\kappa}
#' @seealso \code{\link{dfisher}},\code{\link{rvmises}},\code{\link{rcayley}},\code{\link{rhaar}}


rfisher <- function(n, kappa = 1, nu = NULL) {
  
  if(!is.null(nu))
    kappa <- fisher_kappa(nu)
  
  M <- max(dfisher(seq(-pi, pi, length = 1000), kappa,Haar=F))
  return(rar(n, dfisher, M, kappa = kappa, Haar=F))
}

#' Simulate a data set of size \eqn{n} from the uniform distribution on the sphere
#'
#' The uniform distribution has the density\deqn{C_\mathrm{{F}}(r|\kappa)=\frac{1}{2\pi}1-\cos(r)}.  The is also 
#' know as the Haar measure.
#' 
#' @param n sample size
#' @return a sample of size \eqn{n} from the uniform distribution on the sphere
#' @seealso \code{\link{dhaar}},\code{\link{rfisher}},\code{\link{rvmises}},\code{\link{rcayley}}


rhaar<-function(n){
  return(rar(n, dhaar, 1/pi))
}


#' Generate a vector of angles(r) from the von Mises Circular distribution
#'
#' The circular von Mises-based distribution has the density \deqn{C_\mathrm{M}(r|\kappa)=\frac{1}{2\pi \mathrm{I_0}(\kappa)}e^{\kappa\cos(r)}}.  This function allows the use to
#' simulate \eqn{n} random deviates from \eqn{C_\mathrm{M}(r|\kappa)} given a concentration parameter \eqn{\kappa}.
#'
#' @param n The number of angles desired
#' @param kappa The concentration parameter of the distribution
#' @param nu An alternative to kappa; circular variance
#' @return S3 \code{rvmises} object; a vector of n angles following the von Mises Circular distribution with concentration kappa and mean/mode 0
#' @export
#' @examples
#' r<-rvmises(20,0.01)

rvmises <- function(n, kappa = 1, nu = NULL) {
  
  if(!is.null(nu))
    kappa <- vmises_kappa(nu)
  
  u <- runif(3, 0, 1)
  a <- 1 + sqrt(1 + 4 * kappa^2)
  b <- (a - sqrt(2 * a))/(2 * kappa)
  r <- (1 + b^2)/(2 * b)
  theta <- rep(10, n)
  
  for (i in 1:n) {
    
    while (theta[i] == 10) {
      # Step 1
      u <- runif(3, 0, 1)
      z <- cos(pi * u[1])
      f <- (1 + r * z)/(r + z)
      c <- kappa * (r - f)
      
      # Step 2
      u <- runif(3, 0, 1)
      if ((c * (2 - c) - u[2]) > 0) {
        
        theta[i] = sign(u[3] - 0.5) * acos(f)
        
      } else {
        
        if (log(c/u[2]) + 1 - c < 0) {
          u <- runif(3, 0, 1)
        } else {
          u <- runif(3, 0, 1)
          theta[i] = sign(u[3] - 0.5) * acos(f)
        }
      }
    }
  }
  return(theta)
}


#' Compute the sum of the \eqn{p^{\text{th}}} order distances between Rs and S
#'
#' @param Rs a matrix of rotation observations, one row per observation
#' @param S the individual matrix of interest, usually an estimate of the mean
#' @param p the order of the distances to compute
#' @return the sum of the pth order distance between each sample in Rs and S
#' @export
#' @examples
#' r<-rvmises(20,0.01)
#' Rs<-genR(r)
#' Sp<-mean(Rs)
#' sum_dist(Rs,S=Sp,2)

sum_dist<-function(Rs, S = id.SO3, method='projected', p=1){
  
  UseMethod( "sum_dist" )

}

#' @return \code{NULL}
#'
#' @rdname sum_dist
#' @method sum_dist EA
#' @S3method sum_dist EA

sum_dist.EA <- function(EAs, S = id.EA, method='projected', p=1) {
  
  return(sum(apply(EAs, 1, dist.EA , EA2 = S, method, p)))
  
}

#' @return \code{NULL}
#'
#' @rdname sum_dist
#' @method sum_dist Q4
#' @S3method sum_dist Q4

sum_dist.Q4 <- function(Qs, S = id.Q4, method='projected', p=1) {
  
  return(sum(apply(Qs, 1, dist.Q4 , Q2 = S, method, p)))
  
}


#' @return \code{NULL}
#'
#' @rdname sum_dist
#' @method sum_dist SO3
#' @S3method sum_dist SO3

sum_dist.SO3 <- function(Rs, S = id.SO3, method='projected', p=1) {
  
  return(sum(apply(Rs, 1, dist.SO3 , S = S, method, p)))
  
}


tLogMat <- function(x, S) {
  tra <- log.SO3(t(S) %*% matrix(x, 3, 3))
  return(as.vector(tra))
}


vecNorm <- function(x, S, ...) {
  n <- sqrt(length(x))
  cenX <- x - as.vector(S)
  return(norm(matrix(cenX, n, n), ...))
}

#' Identity in SO(3) space
#' @export
id.SO3 <- as.SO3(diag(c(1,1,1)))

#' Identity in Q4 space
#' @export
id.Q4 <- as.Q4(c(1,0,0,0))

#' Identity in EA space
#' @export
id.EA <- as.EA(c(0,0,0))

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

