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




#' Distance Between Two Rotations
#'
#' This function will calculate the intrinsic (Riemannian) or projected (Euclidean) distance between two rotations.  If only one rotation is specified
#' the other will be set to the identity and the distance between the two is returned.
#'
#' @param x rotation in SO3 representation
#' @param ... Additional arguments
#' @return the distance between x and something else
#' @export

dist<-function(x,...){
  UseMethod("dist")
}

#' @return \code{NULL}
#' 
#' @rdname dist
#' @method dist SO3
#' @S3method dist SO3

dist.SO3 <- function(R1, R2=id.SO3, method='projected' , p=1) {
  
  R1<-as.SO3(matrix(R1,3,3))
  R2<-as.SO3(matrix(R2,3,3))
  
  if(method=='projected'){
    
    so3dist<-norm(R1-R2,type='F')^p
    
  }else if(method=='intrinsic'){
    
    so3dist<-angle(as.SO3(t(R1)%*%R2))^p
    
  }else{
    stop("Incorrect usage of method argument.  Please choose intrinsic or projected.")
  }
  
  return(so3dist)
  
}

#' @return \code{NULL}
#' 
#' @rdname dist
#' @method dist Q4
#' @S3method dist Q4

dist.Q4 <- function(R1, R2=id.Q4 ,method='projected', p=1) {
  Q1 <- R1
  Q2 <- R2
  
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

#' @return \code{NULL}
#' 
#' @rdname dist
#' @method dist EA
#' @S3method dist EA

dist.EA <- function(R1, R2=id.EA ,method='projected', p=1) {
  EA1 <- R1
  EA2 <- R2
  
  R1<-SO3.EA(EA1)
  
  R2<-SO3.EA(EA2)
  
  EAdist<-dist.SO3(R1,R2,method,p)
  
  return(EAdist)
}




#' Find the angle of rotation R
#' 
#' Extract angle from rotation.
#' 
#' @param Rs rotation matrix
#' @return angle of rotation
#' @seealso \code{\link{axis}}
#' @export

angle<-function(Rs){
  UseMethod("angle")
}

#' @return \code{NULL}
#'
#' @rdname angle
#' @method angle SO3
#' @S3method angle SO3

angle.SO3 <- function(Rs){
  ##  trace of a rotation matrix has to be between -1 and 3. If not, this is due
  ## to numerical inconcistencies, that we have to fix here
  tr<-Rs[1]+Rs[5]+Rs[9]
  eps <- 10^-3
  
if (tr > 3+eps) print(sprintf("Warning: trace too large (> 3.0): %f ", tr))
  if (tr < -1-eps) print(sprintf("Warning: trace too small (< -1.0): %f ", tr))
#  stopifnot(tr<3+eps, tr>-1-eps)
  tr <- max(min(3, tr), -1)
  
  return(acos((tr-1)/2))
}

#' @return \code{NULL}
#'
#' @rdname angle
#' @method angle Q4
#' @S3method angle Q4

angle.Q4 <- function(Qs){
  theta<-2*acos(Qs[1])
  return(theta)
}

#' @return \code{NULL}
#'
#' @rdname angle
#' @method angle EA
#' @S3method angle EA

angle.EA<-function(eur){
  
  trR <- cos(eur[1]) * cos(eur[3]) - sin(eur[1]) * sin(eur[3]) * cos(eur[2])
  trR <- trR + (-sin(eur[1]) * sin(eur[3]) + cos(eur[1]) * cos(eur[3]) * cos(eur[2]))
  trR <- trR + (cos(eur[2]))
  theta<-acos((trR-1)/2)
  return(theta)
}


#' Find the axis of rotation R
#' 
#' This function will find the axis of rotation matrix R.  The simple calculation is based on Rodrigues formula
#' and noticing that R - t(R) can be simplified greatly.  
#' @param R 3-by-3 matrix in SO3 
#' @return axis in form of three dimensional vector of length one.
#' @seealso \code{\link{angle}}
#' @export

axis2<-function(R){
  UseMethod("axis2")
}

#' @return \code{NULL}
#'
#' @rdname axis2
#' @method axis2 SO3
#' @S3method axis2 SO3

axis2.SO3<-function(R){
  # based on Rodrigues formula: R - t(R)
  R<-matrix(R,3,3)
  
  X <- R - t(R)
  u <- rev(X[upper.tri(X)])*c(-1,1,-1)
  
  return(u/sqrt(sum(u^2))) # will be trouble, if R is symmetric, i.e. id,  .... 

}

#' @return \code{NULL}
#'
#' @rdname axis2
#' @method axis2 Q4
#' @S3method axis2 Q4

axis2.Q4 <- function(q){
  
  theta<-angle(q)
  
  if(theta==0){
    u <- rep(0,3)
  }else{
    u <- q[2:4]/sin(theta/2)
  }
  return(u)
}

#' @return \code{NULL}
#'
#' @rdname axis2
#' @method axis2 EA
#' @S3method axis2 EA

axis2.EA <- function(eur){
  R<-SO3.EA(eur)
  u<-axis2(R)
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
#' plot(Rs,center=mean(Rs),show_estimates=TRUE,shape=4)

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
    PMed <- as.vector(median.SO3(Rs))
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
    PMed <- as.vector(median.SO3(Rs))
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
    stop("Input has to be of class SO(3).")
  }
  
  theta <- angle.SO3(R)
  
  if (abs(cos(theta)) >= 1) {
    return(diag(0, 3, 3))
  } else {
    
    c <- theta/(2 * sin(theta))
    mlog <- c * (R - t(R))
    return(mlog)
  }
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


#' Compute the sum of the \eqn{p^{\text{th}}} order distances between Rs and S
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

  return(sum(apply(Rs, 1, dist.SO3 , R2 = S, method=method, p=p)))
  
}
#' @return \code{NULL}
#'
#' @rdname sum_dist
#' @method sum_dist EA
#' @S3method sum_dist EA

sum_dist.EA <- function(EAs, S = id.EA, method='projected', p=1) {
  
  return(sum(apply(EAs, 1, dist.EA , EA2 = S, method=method, p=p)))
  
}

#' @return \code{NULL}
#'
#' @rdname sum_dist
#' @method sum_dist Q4
#' @S3method sum_dist Q4

sum_dist.Q4 <- function(Qs, S = id.Q4, method='projected', p=1) {
  
  return(sum(apply(Qs, 1, dist.Q4 , Q2 = S, method=method, p=p)))
  
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




