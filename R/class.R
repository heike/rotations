#' Class of rotation matrices. 
#' 
#' An S4 class that stores a matrix describing a 3 by 3 rotational matrix.
#' @slot R is a 3 by 3 matrix
#' @export
#' @examples
#' a <- new("rotation") ## a is identity
#' b <- new("rotation", R=matrix(runif(9), nrow=3))
#' ## b is not a rotation matrix
setClass("rotation",
  representation(R = "matrix", U="numeric", theta="numeric"),
  prototype(R=diag(c(1,1,1)), U=c(1,0,0), theta=0)
)



#' Method for testing for 3d rotations
#'
#' A 3d rotation matrix is defined as a 3 by 3 dimensional matrix 
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


	R <- P + (id - P) * cos(theta) + eskew(U) * sin(theta)
  new("rotation", R=R, U=U, theta=theta)
}