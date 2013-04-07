library(ggplot2)
require(grid)

# set origin of concentric circles
origin <- matrix(SO3(c(1,-1,0), pi/16),3,3)
  
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


setOrigin <- function(origin = matrix(SO3(c(1,-1,0), pi/8),3,3)) {
  origin <<- origin
  pcircles <- data.frame(as.matrix(circles[,1:3]) %*% origin)
  pcircles
}


# this is the coordinate system and should be fixed, no matter what column of the rotation matrices is shown

base <- ggplot(aes(x=X1, y=X2), data=setOrigin(matrix(SO3(c(1,-1,0), pi/16),3,3))) + 
  coord_equal() + 
  geom_point(aes(alpha=X3), size=0.6, colour="grey65") + 
  scale_alpha(range=c(0,0.8),  guide="none") + 
  theme(panel.background=element_blank(),
       panel.grid.minor=element_blank(),
       panel.grid.major=element_blank(),
       axis.title.x=element_blank(),
       axis.title.y=element_blank(),
       axis.text.x=element_blank(),
       axis.text.y=element_blank(),
       axis.ticks=element_blank(), 
       plot.margin = unit(rep(0, 4), "lines"))


roteye <- function(origin, center, column=1) {
  R <- list(matrix(SO3(c(0,1,0), pi/2),3,3), matrix(SO3(c(1,0,0), -pi/2),3,3), diag(c(1,1,1)))[[column]]
  rot <- center %*% R %*% origin 
}


#' Project rotation data onto sphere
#' 
#' Projection of rotation matrices onto sphere with given center.
#'
#' @param data data frame of rotation matrices in 3 x 3 matrix representation
#' @param center point about which to center the observations
#' @param column integer 1 to 3 indicating which column to display
#' @return  data frame with columns X, Y, Z standing for the respective coordinates in 3d space
#' @export
#' 
pointsXYZ <- function(data, center, column=1) {
  rot <- roteye(origin, center, column)
  idx <- list(1:3,4:6, 7:9)[[column]]
  data <- as.matrix(data[,idx])
  
  psample1 <- data.frame(data %*% rot)
  names(psample1) <- c("X","Y","Z")
  
  #  psample1 <- data.frame(psample1, data)
  #  psample1 <- psample1[order(psample1$Z, decreasing=FALSE),]
  psample1  
}


#' Visualizing random rotations.
#'
#' This function produces a three-dimensional globe onto which  one of the  columns of the provided sample of rotations is drawn.  The data are centered around a provided
#' matrix and the user can choose to display this center or not.  Based on \code{ggplot2} package by \cite{wickham09}.
#'
#' @param x n rotations in SO3 format
#' @param center point about which to center the observations
#' @param col integer 1 to 3 indicating which column to display
#' @param toRange show only part of the globe that is in range of the data?
#' @param show_estimates character vector to specify  which of the four estimates of the principal direction to show. Possibilities are
#'     "all", "proj.mean", "proj.median", "riem.mean", "riem.median"
#' @param ... parameters passed onto the points layer
#' @return  a ggplot2 object with the data displayed on spherical grid
#' @cite wickham09
#' @export
#' @examples
#' r<-rvmises(200,1.0)
#' Rs<-genR(r)
#' plot(Rs,center=mean(Rs),show_estimates=NULL,shape=4)
#' # Z is computed internally and contains information on depth
#' plot(Rs,center=mean(Rs),show_estimates=c("proj.mean", "riem.mean")) + aes(size=Z, alpha=Z) + scale_size(limits=c(-1,1), range=c(0.5,2.5))

plot.SO3 <- function(x, center, col=1, toRange=FALSE, show_estimates=NULL,  ...) {
  Rs <- as.SO3(x)
  xlimits <- c(-1,1)
  ylimits <- c(-1,1)
  
  X <- Y <- Est <- NULL
  proj2d <- pointsXYZ(Rs, center=center, column=col)
  if(toRange) {
    xlimits <- range(proj2d$X)
    ylimits <- range(proj2d$Y)
    xbar <- mean(xlimits)
    xlimits <- xbar + 1.1*(xlimits-xbar)
    ybar <- mean(ylimits)
    ylimits <- ybar + 1.1*(ylimits-ybar)
  }
  estimates <- NULL
  if (!is.null(show_estimates)) {
    ShatP <- StildeP <- ShatG <- StildeG <- NA
    if(show_estimates%in%c('all','All')) show_estimates<-c("proj.mean","proj.median","riem.mean","riem.median")
    if (length(grep("proj.mean", show_estimates)) > 0) ShatP<-mean(Rs, type="projected")
    if (length(grep("proj.median", show_estimates)) >0)    StildeP<-median(Rs, type="projected")
    if (length(grep("riem.mean", show_estimates)) > 0)    ShatG<-mean(Rs, type="intrinsic")
    if (length(grep("riem.median", show_estimates)) > 0)    StildeG<-median(Rs, type="intrinsic")
    
    Shats<-data.frame(rbind(as.vector(ShatP),as.vector(StildeP),as.vector(ShatG),as.vector(StildeG)),Est=1:4)
    Shats$Est <- factor(Shats$Est)
    labels <- c(expression(hat(S)[E]), expression(tilde(S)[E]), expression(hat(S)[R]), expression(tilde(S)[R]))
    levels(Shats$Est) <- labels
    Shats <- na.omit(Shats)
    
    estimates <- list(geom_point(aes(x=X, y=Y, colour=Est),size=3.5, data=data.frame(pointsXYZ(Shats, center=center, column=col), Shats)),
                      scale_colour_brewer("Estimates", palette="Paired", labels=labels))
  }
  base + geom_point(aes(x=X, y=Y), data=proj2d, ...) + 
    estimates +
    xlim(xlimits) + ylim(ylimits) 
}