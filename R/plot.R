library(ggplot2)
require(grid)
# set origin of concentric circles
origin <- SO3(c(1,-1,0), pi/16)
  
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


setOrigin <- function(origin = SO3(c(1,-1,0), pi/8)) {
  origin <<- origin
  pcircles <- data.frame(as.matrix(circles[,1:3]) %*% origin)
  pcircles
}


# this is the coordinate system and should be fixed, no matter what column of the rotation matrices is shown

base <- ggplot(aes(x=X1, y=X2), data=setOrigin(SO3(c(1,-1,0), pi/16))) + 
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
  R <- list(SO3(c(0,1,0), pi/2), SO3(c(1,0,0), -pi/2), diag(c(1,1,1)))[[column]]
  rot <- center %*% R %*% origin 
}

pointsXY <- function(data, center, column=1) {
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
#' @param Rs the sample of n random rotations
#' @param center point about which to center the observations
#' @param column integer 1 to 3 indicating which column to display
#' @param toRange show only part of the globe that is in range of the data?
#' @param show_estimates rather to display the four estimates of the principal direction or not
#' @return  a ggplot2 object with the data dispalyed on spherical grid
#' @cite wickham09
#' @export
#' @examples
#' r<-rvmises(20,1.0)
#' Rs<-genR(r)
#' eyeBall(Rs,center=mean(Rs),show_estimates=TRUE,shape=4)
#' @export
plot.SO3 <- function(Rs, center, col=1, toRange=FALSE, show_estimates=FALSE) {
  Rs <- as.SO3(Rs)
  xlimits <- c(-1,1)
  ylimits <- c(-1,1)
  
  proj2d <- pointsXY(Rs, center=center, column=col)
  if(toRange) {
    xlimits <- range(proj2d$X)
    ylimits <- range(proj2d$Y)
    xbar <- mean(xlimits)
    xlimits <- xbar + 1.1*(xlimits-xbar)
    ybar <- mean(ylimits)
    ylimits <- ybar + 1.1*(ylimits-ybar)
  }
  estimates <- NULL
  if (show_estimates) {
    ShatP<-mean(Rs, type="projected")
    StildeP<-median.SO3(Rs, type="projected")
    ShatG<-mean(Rs, type="intrinsic")
    StildeG<-median.SO3(Rs, type="intrinsic")
    
    Shats<-data.frame(rbind(as.vector(ShatP),as.vector(StildeP),as.vector(ShatG),as.vector(StildeG)),Est=1:4)
    Shats$Est <- factor(Shats$Est)
    labels <- c(expression(hat(S)[P]), expression(tilde(S)[P]), expression(hat(S)[G]), expression(tilde(S)[G]))
    levels(Shats$Est) <- labels

    estimates <- list(geom_point(aes(x=X, y=Y, colour=Est), size=3, data=data.frame(pointsXY(Shats, center=center, column=col), Shats)),
    scale_colour_brewer("Estimates", palette="Paired", labels=labels))
  }
  base + geom_point(aes(x=X, y=Y), size=2, data=proj2d, colour="grey20") + 
    estimates+
    xlim(xlimits) + ylim(ylimits)
}