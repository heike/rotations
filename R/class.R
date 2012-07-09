check_rotation <- function(object) {
  if (is.null(object@R))
    object@R <- angle_axis(object@theta, object@U)

  if (is.rotation(object)) return(TRUE)
  else return("R is not a rotation matrix. Specify angle and axis instead")
}

#' Class of rotation matrices. 
#' 
#' An S4 class that stores a matrix describing a 3 by 3 rotational matrix.
#' @slot R is a 3 by 3 matrix
#' @slot U is a 3 dimensional vector describing the fix axis of rotation R
#' @slot thetais the angle by which rotation R rotates around U.
#' @export
#' @examples
#' a <- new("rotation") ## a is identity
#' b <- new("rotation", R=matrix(runif(9), nrow=3))
#' ## b is not a rotation matrix
setClass("rotation",
  representation(R = "matrix", U="numeric", theta="numeric"),
  prototype(R=diag(c(1,1,1)), U=NA_real_, theta=NA_real_),
  validity=check_rotation
)



#' Method for testing for 3d rotations
#'
#' A 3 dimensional rotation is defined as a 3 by 3 dimensional matrix R with the following properties:
#' \eqn{R^t R = I_{3 \times 3}}.
#'
#' @export
#' @docType methods
#' @rdname is.rotation-methods
setGeneric("is.rotation", function(object, tol=1e-6) {
  if (!(all(dim(object@R) == c(3,3)))) return(FALSE)
  
  all(abs(object@R %*% t(object@R) - diag(c(1,1,1))) < tol)
})



#' Method for creating rotations using the angle axis representation
#'
#' Angle-axis representation based on the Rodrigues formula.
#'
#' @export
#' @param U three-dimensional vector describing the fix axis of the rotation
#' @param theta angle between -pi and pi
angle_axis <- function(U, theta) {
	U <- U/sqrt(sum(U^2))
	P <- U %*% t(U)

	id <- matrix(0, length(U), length(U))
	diag(id) <- 1


	R <- P + (id - P) * cos(theta) + Phi(U) * sin(theta)
  new("rotation", R=R, U=U, theta=theta)
}

#' Skew-symmetric matrix corresponding to rotation R
#' 
#' The space of skew-symmetric matrices, i.e. all square matrices A with \deqn{A^t = -A}
#' is the tangent space of the special orthogonal group \deqn{SO(3)}, that describes all three-dimensional
#' rotations.
#' @export
#' @param U three-dimensional vector of real numbers describing the fix axis of rotation R.
Phi <- function(U) {
  U <- U/sqrt(sum(U^2))
  u <- U[1]
  v <- U[2]
  w <- U[3]
  
  res <- matrix((-1) *c(0,-w,v,w,0,-u,-v,u,0), ncol=3)
  return(res)
}
